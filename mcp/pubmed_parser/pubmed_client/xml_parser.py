"""
XML and HTML parsing functions for PubMed/PMC articles.
"""

import re
import xml.etree.ElementTree as ET


# ── 标准三字母 -> 单字母氨基酸映射 ──
AA_3TO1 = {
    "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C",
    "Gln": "Q", "Glu": "E", "Gly": "G", "His": "H", "Ile": "I",
    "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F", "Pro": "P",
    "Ser": "S", "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V",
    "Ter": "*", "ter": "*", "Xaa": "X", "X": "*",
    "*": "*",
}

AA_1TO3 = {v: k for k, v in AA_3TO1.items() if v != "X"}
AA_1TO3["*"] = "Ter"

# ── 已知变异历史命名 → HGVS 映射 ──
# 格式: { 基因符号: { 历史名称: {"cdna": "c.xxx", "protein": "p.xxx"} } }
KNOWN_VARIANT_NAMES = {
    "G6PD": {
        "Tsukui": {"cdna": "c.565_567del", "protein": "p.Ser189del"},
        "G6PD Tsukui": {"cdna": "c.565_567del", "protein": "p.Ser189del"},
    },
    # 可继续扩展其他基因的历史命名
}


def _empty_result():
    return {
        "PMID": "", "标题": "", "摘要": "", "全文": "",
        "全文来源": "abstract", "作者": [], "期刊": "",
        "发表年份": "", "MeSH术语": [], "一句话概括": "",
        "基因": "", "cDNA变异": "", "蛋白变异": "",
        "变异类型": "", "致病性": "", "合子状态": "",
        "患者详情": [], "临床表型": "",
        "遗传模式": "", "功能验证": "",
        "功能验证详情": [], "变异提及": False,
        "相关句子": [], "匹配关键词": [],
        "共存变异": [], "反式确认": False,
        # v7 新增;正反式/相位字段
        "顺式确认": False,
        "相位状态": "",           # confirmed_in_trans / confirmed_in_cis / presumed_in_trans / presumed_in_cis / phase_not_determined / not_assessed / not_applicable
        "相位置信度": "",         # confirmed / presumed / unknown / not_applicable
        "亲本检测": False,
        "母源变异": None,
        "父源变异": None,
        "相位证据": [],           # list of {type, label, sentence}
        "患者数量": 0, "变异特征": {},
        "表格数量": 0, "表格摘要": [],
        "总结段落": "", "临床详情": {},
    }


def parse_ncbi_xml(xml_text):
    result = _empty_result()
    if not xml_text:
        return None
    try:
        root = ET.fromstring(xml_text)
    except ET.ParseError:
        return None

    pmid_elem = root.find(".//PMID")
    if pmid_elem is not None:
        result["PMID"] = pmid_elem.text.strip()

    title_elem = root.find(".//ArticleTitle")
    if title_elem is not None:
        result["标题"] = "".join(title_elem.itertext()).strip()

    abstract_parts = []
    for elem in root.findall(".//AbstractText"):
        label = elem.get("Label", "")
        text = "".join(elem.itertext()).strip()
        abstract_parts.append(f"[{label}] {text}" if label else text)
    result["摘要"] = " ".join(abstract_parts)
    result["全文"] = result["摘要"]
    result["全文来源"] = "abstract"

    for author_elem in root.findall(".//Author"):
        ln = author_elem.find("LastName")
        fn = author_elem.find("ForeName")
        parts = []
        if ln is not None:
            parts.append(ln.text or "")
        if fn is not None:
            parts.append(fn.text or "")
        if parts:
            result["作者"].append(" ".join(parts))

    journal_elem = root.find(".//Journal/Title")
    if journal_elem is not None:
        result["期刊"] = journal_elem.text.strip()
    else:
        abbrev = root.find(".//Journal/ISOAbbreviation")
        if abbrev is not None:
            result["期刊"] = abbrev.text.strip()

    pub_date = root.find(".//PubDate/Year")
    if pub_date is not None:
        result["发表年份"] = pub_date.text.strip()
    else:
        medline = root.find(".//PubDate/MedlineDate")
        if medline is not None and medline.text:
            m = re.search(r'(\d{4})', medline.text)
            if m:
                result["发表年份"] = m.group(1)

    return result


def parse_pmc_html(html_text, result):
    """从 PMC HTML 页面提取正文文本和表格数据."""
    if not html_text or not result:
        return result

    # HTML 实体解码
    html_text = html_text.replace("&gt;", ">").replace("&lt;", "<").replace("&amp;", "&")
    html_text = html_text.replace("&quot;", "\"").replace("&#x0002C;", ",")
    html_text = html_text.replace("&nbsp;", " ")

    text_parts = []

    # 1. 提取 abstract
    abs_match = re.search(r'<div[^>]*class="[^"]*abstract[^"]*"[^>]*>(.*?)</div>', html_text, re.DOTALL | re.IGNORECASE)
    if not abs_match:
        abs_match = re.search(r'<section[^>]*id="[^"]*abstract[^"]*"[^>]*>(.*?)</section>', html_text, re.DOTALL | re.IGNORECASE)
    if abs_match:
        abs_text = re.sub(r'<[^>]+>', ' ', abs_match.group(1))
        abs_text = re.sub(r'\s+', ' ', abs_text).strip()
        if abs_text:
            text_parts.append(abs_text)

    # 2. 从 <section class="body main-article-body"> 提取正文(PMC 标准结构)
    body_match = re.search(
        r'<section[^>]*class="[^"]*\bbody\b[^"]*"[^>]*>(.*?)</section>',
        html_text, re.DOTALL | re.IGNORECASE
    )
    if body_match:
        # 排除表格部分(表格单独提取)
        body_html = re.sub(r'<table[^>]*>.*?</table>', '', body_match.group(1), flags=re.DOTALL | re.IGNORECASE)
        # 提取所有 <p> 标签文本
        paras = re.findall(r'<p[^>]*>(.*?)</p>', body_html, re.DOTALL)
        for p in paras:
            p_text = re.sub(r'<[^>]+>', ' ', p)
            p_text = re.sub(r'\s+', ' ', p_text).strip()
            if len(p_text) > 50:
                text_parts.append(p_text)

    # 3. 通用 fallback;从所有 <p> 标签提取(排除表头/页脚杂讯)
    if not text_parts:
        all_paras = re.findall(r'<p[^>]*>(.*?)</p>', html_text, re.DOTALL)
        for p in all_paras:
            p_text = re.sub(r'<[^>]+>', ' ', p)
            p_text = re.sub(r'\s+', ' ', p_text).strip()
            # 排除 PMC 界面杂讯
            if len(p_text) > 80 and not re.match(
                r'(An official website|The \.gov|Federal government|The site is secure|'
                r'Access keys|NCBI Homepage|MyNCBI|PubMed|PMC|Follow|Share|Connect|'
                r'Disclaimer|Copyright|FOIA|Privacy|NLM|National Library|NIH|HHS|'
                r'Vulnerability Disclosure|Accessibility|Careers|Nondiscrimination)',
                p_text, re.IGNORECASE
            ):
                text_parts.append(p_text)

    if text_parts:
        full_text = " ".join(text_parts)
        if len(full_text) > len(result.get("全文", "")):
            result["全文"] = full_text
            result["全文来源"] = "pmc_html"

    # 提取表格数据
    tables = _extract_tables_from_html(html_text)
    if tables:
        result["tables"] = tables

    return result


def _extract_tables_from_html(html_text):
    """从 HTML 文本提取表格(正则方式,零依赖)."""
    tables = []
    # 找所有 table 元素
    table_matches = list(re.finditer(r'<table[^>]*>(.*?)</table>', html_text, re.DOTALL | re.IGNORECASE))
    for table_match in table_matches:
        table_content = table_match.group(1)
        # 查找 caption(可能在 table 之前或内部)
        caption = ""
        cap_match = re.search(r'<caption[^>]*>(.*?)</caption>', table_content, re.DOTALL | re.IGNORECASE)
        if cap_match:
            caption = re.sub(r'<[^>]+>', ' ', cap_match.group(1)).strip()
            caption = re.sub(r'\s+', ' ', caption)

        rows = []
        tr_matches = re.finditer(r'<tr[^>]*>(.*?)</tr>', table_content, re.DOTALL | re.IGNORECASE)
        for tr in tr_matches:
            tr_content = tr.group(1)
            cells = []
            for cell_match in re.finditer(r'<(?:td|th)[^>]*>(.*?)</(?:td|th)>', tr_content, re.DOTALL | re.IGNORECASE):
                cell_text = re.sub(r'<[^>]+>', ' ', cell_match.group(1))
                cell_text = re.sub(r'\s+', ' ', cell_text).strip()
                if re.match(r'<(?:td|th)[^>]*>', tr_content[:50]):
                    # 标记 header
                    cell_text = f"[H]{cell_text}"
                cells.append(cell_text)
            if cells:
                rows.append(cells)

        if rows:
            tables.append({
                "id": "",
                "caption": caption,
                "rows": rows,
                "raw_text": "\n".join(
                    "\t".join(c for c in row) for row in rows
                ),
            })
    return tables


def parse_pmc_europe_xml(xml_text, result):
    if not xml_text or not result:
        return result

    # 检测出版商限制(如 ASN 出版社会在 XML 中标注限制)
    if "does not allow downloading of the full text in XML form" in xml_text:
        try:
            root = ET.fromstring(xml_text)
        except ET.ParseError:
            return result
        # 提取 front 部分的文本(标题,摘要等)
        front_texts = []
        for elem in root.findall(".//front//*"):
            if elem.text and len(elem.text.strip()) > 20:
                front_texts.append(elem.text.strip())
        if front_texts:
            result["全文"] = " ".join(front_texts)
            result["全文来源"] = "pmc_restricted"
        return result

    try:
        root = ET.fromstring(xml_text)
    except ET.ParseError:
        return result

    # 提取全文文本(优先 body,同时提取摘要)
    text_parts = []

    # 1. 先提取 PMC XML 中的摘要(front 部分)
    for abs_elem in root.findall(".//abstract/abstract-text"):
        txt = " ".join(abs_elem.itertext()).strip()
        if txt:
            text_parts.append(txt)
    if not text_parts:
        for abs_elem in root.findall(".//AbstractText"):
            txt = "".join(abs_elem.itertext()).strip()
            if txt:
                text_parts.append(txt)

    # 2. 提取 body 文本
    body_texts = []
    body = root.find(".//body")
    if body is not None:
        for elem in body.iter():
            if elem.text:
                body_texts.append(elem.text.strip())
            if elem.tail:
                body_texts.append(elem.tail.strip())

    # 3. 如果 body 不存在或为空,尝试从 <sec> 元素提取(非JATS标准结构)
    if not body_texts:
        for sec in root.findall(".//sec"):
            for elem in sec.iter():
                if elem.text and len(elem.text.strip()) > 10:
                    body_texts.append(elem.text.strip())
                if elem.tail and len(elem.tail.strip()) > 10:
                    body_texts.append(elem.tail.strip())
    if not body_texts:
        for elem in root.iter():
            if elem.tag not in ("article", "front", "body", "back"):
                if elem.text:
                    body_texts.append(elem.text.strip())
                if elem.tail:
                    body_texts.append(elem.tail.strip())

    # 合并;摘要 + body
    if body_texts:
        text_parts.extend(body_texts)

    # 如果 PMC 完全没有文本,保留已有的摘要
    if not text_parts:
        return result

    full_text = " ".join(filter(None, text_parts))
    if full_text:
        result["全文"] = full_text
        result["全文来源"] = "pmc_fulltext"

    # 提取表格数据
    tables = _extract_tables_from_xml(root)
    if tables:
        result["tables"] = tables

    return result


def _extract_tables_from_xml(root):
    """从 XML 全文提取表格数据(行列表).支持 table-wrap (JATS) 和 table 元素."""
    tables = []

    def _process_table_element(table_elem, caption=""):
        """处理单个表格元素(可能是 <table> 或 <table-wrap> 内的 <table>)"""
        rows = []
        for tr in table_elem.findall(".//tr"):
            cells = []
            for cell in tr.findall("th"):
                txt = " ".join(t.strip() for t in cell.itertext() if t.strip())
                cells.append(f"[H]{txt}")
            for cell in tr.findall("td"):
                txt = " ".join(t.strip() for t in cell.itertext() if t.strip())
                cells.append(txt)
            if cells:
                rows.append(cells)
        return rows

    # JATS 格式: <table-wrap> 包裹 <table>,caption 在 table-wrap 层
    for tw in root.findall(".//table-wrap"):
        tw_id = tw.get("id", "")
        caption = ""
        caption_el = tw.find("caption")
        if caption_el is not None:
            caption = " ".join(
                t.strip() for t in caption_el.itertext() if t.strip()
            )
        # 在 table-wrap 内找 <table>
        table_el = tw.find("table")
        if table_el is not None:
            rows = _process_table_element(table_el, caption)
        else:
            # table-wrap 内可能直接有 <tr> (某些格式)
            rows = _process_table_element(tw, caption)
        if rows:
            tables.append({
                "id": tw_id,
                "caption": caption,
                "rows": rows,
                "raw_text": "\n".join(
                    "\t".join(c for c in row) for row in rows
                ),
            })

    # 直接的 <table> 元素(不在 table-wrap 内的,如 HTML 格式)
    processed_tables = set()
    for tw in root.findall(".//table-wrap"):
        table_el = tw.find("table")
        if table_el is not None:
            processed_tables.add(table_el)

    for table_elem in root.findall(".//table"):
        if table_elem in processed_tables:
            continue  # 已通过 table-wrap 处理
        caption_el = table_elem.find(".//caption")
        caption = ""
        if caption_el is not None:
            caption = " ".join(
                t.strip() for t in caption_el.itertext() if t.strip()
            )
        rows = _process_table_element(table_elem, caption)
        if rows:
            tables.append({
                "id": table_elem.get("id", ""),
                "caption": caption,
                "rows": rows,
                "raw_text": "\n".join(
                    "\t".join(c for c in row) for row in rows
                ),
            })

    return tables


# ── 句子拆分 ──

def split_sentences(text):
    """split sentences, protect c.NNN format from truncation."""
    protected = re.sub(r'\bc\.(?=\d)', 'c<DOT>', text)
    raw_sentences = re.split(r'(?<=[.!?])\s+', protected)
    return [s.replace('c<DOT>', 'c.').strip() for s in raw_sentences if len(s.strip()) >= 15]