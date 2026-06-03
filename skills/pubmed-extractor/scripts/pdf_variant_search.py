#!/usr/bin/env python3
"""
离线 PDF 变异搜索器 — 确定性、可重复的 PDF 全文变异搜索。
仅在本地 PDF 文件中搜索目标变异，不联网。

用法:
  python pdf_variant_search.py \
      --pdf-dir D:\claude_code\project1\sci \
      --gene ABCA4 \
      --variant "c.763C>T (p.Arg255Cys)" \
      [--transcript NM_000350.3] \
      [--output results.json] \
      [--excel-dir D:\claude_code\project1\文献提取结果]
"""

import argparse
import json
import os
import re
import sys
import csv

# Windows UTF-8 兼容
if sys.platform == 'win32':
    sys.stdout.reconfigure(encoding='utf-8')
    sys.stderr.reconfigure(encoding='utf-8')

try:
    import pdfplumber
except ImportError:
    print("错误: 需要安装 pdfplumber，请运行: pip install pdfplumber")
    sys.exit(1)


# ==================== 变异关键词生成（确定性） ====================

def generate_variant_keywords(variant_str):
    """
    从输入的变异字符串中生成所有搜索关键词变体。
    输入: "c.763C>T (p.Arg255Cys)" 或 "c.1761-2A>G (p.?), intron 12"
    输出: {"exact_cdna": [...], "fuzzy_cdna": [...], "protein": [...]}
    """
    keywords = {"exact_cdna": [], "fuzzy_cdna": [], "protein": []}

    # 解析 cDNA 部分 — 支持多种格式
    # 标准错义: c.763C>T
    # 内含子剪接: c.1761-2A>G, c.3051-1G>A, c.5584+5G>A
    # 缺失: c.101_106delCTTTAT
    # 插入: c.2063_2064insA
    # 重复: c.XXXXdup

    # v8: 通用 cDNA 匹配 — 匹配 c. 后跟数字、字母、符号的组合
    cdna_generic = re.search(r'c\.([\d\w_\-+>delinsup]+)', variant_str.replace(" ", ""))
    if cdna_generic:
        cdna_full = f"c.{cdna_generic.group(1)}"
        keywords["exact_cdna"].append(cdna_full)
        # 提取纯数字部分用于模糊匹配（仅保留 >=3 位的数字，避免 1-2 位数字造成大量假阳性）
        nums = re.findall(r'(\d+)', cdna_generic.group(1))
        for n in nums:
            if n not in keywords["fuzzy_cdna"] and len(n) >= 3:
                keywords["fuzzy_cdna"].append(n)
        # 去除 c. 前缀
        bare = cdna_generic.group(1)
        keywords["fuzzy_cdna"].append(bare)
        keywords["fuzzy_cdna"].append(f"c{bare}")

        # 内含子变异简化形式: 1761-2A>G → 也添加 "2A>G" 形式的变体
        if '-' in bare or '+' in bare:
            # 提取数字后的部分: "1761-2A>G" → "2A>G"
            suffix_match = re.search(r'[-\+](\d+[A-ZTCG]>[A-ZTCG])', bare)
            if suffix_match:
                suffix = suffix_match.group(1)
                keywords["fuzzy_cdna"].append(suffix)

    # 解析蛋白部分 — 支持 p.? 未知蛋白变化
    # 标准格式: p.Arg255Cys, p.R255C
    prot_match = re.search(r'p\.\(?([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})(?:Ter|fs|\*)?\)?', variant_str)
    if not prot_match:
        # 单字母格式: p.R255C, p.R255*
        prot_match = re.search(r'p\.\(?([A-Z\*])(\d+)([A-Z\*\?])\)?', variant_str)
    if not prot_match:
        # p.? 未知蛋白变化
        prot_match_q = re.search(r'p\.\s*\?', variant_str)
        if prot_match_q:
            keywords["protein"].extend(["p.?", "?"])

    if prot_match and len(prot_match.groups()) >= 3 and prot_match.group(3) != '?':
        groups = prot_match.groups()
        if len(groups[0]) == 3:
            # 三字母格式
            aa_ref, aa_pos, aa_alt = groups[0], groups[1], groups[2]
            keywords["protein"].extend([
                f"p.{aa_ref}{aa_pos}{aa_alt}",
                f"{aa_ref}{aa_pos}{aa_alt}",
            ])
            # 单字母形式
            aa1_map = {
                'Ala': 'A', 'Arg': 'R', 'Asn': 'N', 'Asp': 'D', 'Cys': 'C',
                'Gln': 'Q', 'Glu': 'E', 'Gly': 'G', 'His': 'H', 'Ile': 'I',
                'Leu': 'L', 'Lys': 'K', 'Met': 'M', 'Phe': 'F', 'Pro': 'P',
                'Ser': 'S', 'Thr': 'T', 'Trp': 'W', 'Tyr': 'Y', 'Val': 'V',
                'Ter': '*', 'Stop': '*'
            }
            ref_1 = aa1_map.get(aa_ref, aa_ref)
            alt_1 = aa1_map.get(aa_alt, aa_alt)
            keywords["protein"].extend([
                f"p.{ref_1}{aa_pos}{alt_1}",
                f"{ref_1}{aa_pos}{alt_1}",
            ])
            # 带空格变体
            keywords["protein"].append(f"{ref_1} {aa_pos} {alt_1}")

            # 终止密码子变体
            if aa_alt in ('Ter', 'Stop', '*'):
                keywords["protein"].extend([f"{ref_1}{aa_pos}*", f"{ref_1}{aa_pos}Ter", f"{ref_1}{aa_pos}X"])
        else:
            # 单字母格式
            ref_1, pos, alt_1 = groups[0], groups[1], groups[2]
            aa1_map_rev = {'A': 'Ala', 'R': 'Arg', 'N': 'Asn', 'D': 'Asp', 'C': 'Cys',
                          'Q': 'Gln', 'E': 'Glu', 'G': 'Gly', 'H': 'His', 'I': 'Ile',
                          'L': 'Leu', 'K': 'Lys', 'M': 'Met', 'F': 'Phe', 'P': 'Pro',
                          'S': 'Ser', 'T': 'Thr', 'W': 'Trp', 'Y': 'Tyr', 'V': 'Val',
                          '*': 'Ter'}
            ref_3 = aa1_map_rev.get(ref_1, ref_1)
            alt_3 = aa1_map_rev.get(alt_1, alt_1)
            keywords["protein"].extend([
                f"p.{ref_1}{pos}{alt_1}",
                f"{ref_1}{pos}{alt_1}",
                f"p.{ref_3}{pos}{alt_3}",
                f"{ref_3}{pos}{alt_3}",
            ])

    # 去重并保持顺序
    for k in keywords:
        seen = set()
        unique = []
        for item in keywords[k]:
            if item not in seen:
                seen.add(item)
                unique.append(item)
        keywords[k] = unique

    return keywords


# ==================== PDF 文本提取 ====================

def extract_pdf_text(pdf_path):
    """使用 pdfplumber 提取 PDF 全文文本。返回 (文本, 表格列表)。"""
    text_parts = []
    tables = []
    try:
        with pdfplumber.open(pdf_path) as pdf:
            for page in pdf.pages:
                t = page.extract_text()
                if t:
                    text_parts.append(t)
                # 提取表格
                page_tables = page.extract_tables()
                for tbl in page_tables:
                    if tbl:
                        tables.append(tbl)
        return '\n'.join(text_parts), tables
    except Exception as e:
        print(f"  警告: PDF 提取失败 {pdf_path}: {e}")
        return "", []


# ==================== 变异搜索 ====================

def search_variant_in_text(text, keywords, gene):
    """
    在文本中搜索目标变异关键词。
    返回: (是否提及, 匹配关键词列表, 相关句子列表)
    """
    all_keywords = []
    for cat in ["exact_cdna", "fuzzy_cdna", "protein"]:
        all_keywords.extend(keywords.get(cat, []))

    matched_keywords = []
    matched_sentences = []

    for kw in all_keywords:
        # 使用大小写不敏感搜索
        pattern = re.compile(re.escape(kw), re.IGNORECASE)
        for m in pattern.finditer(text):
            if kw not in matched_keywords:
                matched_keywords.append(kw)
            # 提取上下文 (前后 200 字符)
            start = max(0, m.start() - 200)
            end = min(len(text), m.end() + 200)
            ctx = text[start:end].replace('\n', ' ').replace('\r', ' ').strip()
            # 尝试提取完整句子
            sentence = extract_sentence(text, m.start(), m.end())
            matched_sentences.append(sentence)

    # 去重句子
    seen = set()
    unique_sentences = []
    for s in matched_sentences:
        if s not in seen:
            seen.add(s)
            unique_sentences.append(s)

    return len(matched_keywords) > 0, matched_keywords, unique_sentences[:20]


def extract_sentence(text, start, end):
    """提取包含匹配位置的完整句子。"""
    # 向前找句首
    sent_start = start
    for i in range(start, max(0, start - 500), -1):
        if text[i] in '.!\n' and i < start - 5:
            sent_start = i + 1
            break
    # 向后找句尾
    sent_end = end
    for i in range(end, min(len(text), end + 500)):
        if text[i] in '.!\n' and i > end + 5:
            sent_end = i + 1
            break
    return text[max(0, sent_start):sent_end].replace('\n', ' ').replace('\r', ' ').strip()


# ==================== 表格搜索 ====================

def search_variant_in_tables(tables, keywords, gene):
    """在表格中搜索目标变异关键词。返回匹配的表格行列表。"""
    all_keywords = []
    for cat in ["exact_cdna", "fuzzy_cdna", "protein"]:
        all_keywords.extend(keywords.get(cat, []))

    matched_rows = []
    for tbl_idx, table in enumerate(tables):
        for row_idx, row in enumerate(table):
            row_text = ' '.join(str(cell) for cell in row if cell)
            for kw in all_keywords:
                if re.search(re.escape(kw), row_text, re.IGNORECASE):
                    matched_rows.append({
                        "table_index": tbl_idx,
                        "row_index": row_idx,
                        "row_content": row,
                        "matched_keyword": kw,
                    })
                    break
    return matched_rows


# ==================== 表型信息提取 ====================

PHENOTYPE_KEYWORDS = [
    "Stargardt", "STGD", "cone-rod", "cone rod", "CRD",
    "retinitis pigmentosa", "RP", "macular dystrophy", "MD",
    "retinal dystrophy", "retinal degeneration", "IRD",
    "orofacial cleft", "cleft lip", "cleft palate", "CLP", "NSCLP",
    "fundus flavimaculatus", "FFM", "AMD",
]

ZYGOSITY_KEYWORDS = {
    "homozygous": "纯合",
    "heterozygous": "杂合",
    "compound heterozygous": "复合杂合",
    "biallelic": "双等位基因",
    "monoallelic": "单等位基因",
    "hemizygous": "半合子",
}

PATHOGENICITY_KEYWORDS = {
    "pathogenic": "致病",
    "likely pathogenic": "可能致病",
    "disease-causing": "致病",
    "disease causing": "致病",
    "benign": "良性",
    "likely benign": "可能良性",
    "VUS": "意义未明",
    "uncertain significance": "意义未明",
    "polymorphism": "多态性",
    "probably damaging": "可能有害",
    "possibly damaging": "可能有害",
    "tolerated": "耐受",
    "deleterious": "有害",
}

PHASE_KEYWORDS = {
    "in trans": ("confirmed_in_trans", "确认反式"),
    "trans configuration": ("confirmed_in_trans", "确认反式"),
    "on opposite alleles": ("confirmed_in_trans", "确认反式"),
    "on different alleles": ("confirmed_in_trans", "确认反式"),
    "in cis": ("confirmed_in_cis", "确认顺式"),
    "cis configuration": ("confirmed_in_cis", "确认顺式"),
    "on the same allele": ("confirmed_in_cis", "确认顺式"),
    "complex allele": ("confirmed_in_cis", "确认顺式"),
    "compound heterozygous": ("presumed_in_trans", "推定反式"),
    "phase not determined": ("phase_not_determined", "相位未确定"),
    "phase unknown": ("phase_not_determined", "相位未确定"),
}


def extract_phenotype_info(text, sentences):
    """从文本中提取表型、合子状态、致病性信息。"""
    text_lower = text.lower()

    # 表型
    phenotypes = []
    for kw in PHENOTYPE_KEYWORDS:
        if kw.lower() in text_lower:
            if kw not in phenotypes:
                phenotypes.append(kw)

    # 合子状态
    zygosity = []
    for eng, chn in ZYGOSITY_KEYWORDS.items():
        if eng in text_lower:
            zygosity.append(chn)

    # 致病性（从匹配的句子中查找）
    pathogenicity = []
    search_text = ' '.join(sentences).lower() if sentences else text_lower
    for eng, chn in PATHOGENICITY_KEYWORDS.items():
        if eng.lower() in search_text:
            pathogenicity.append(chn)

    # 相位信息
    phase_state = "not_assessed"
    phase_confidence = ""
    for eng, (state, label) in PHASE_KEYWORDS.items():
        if eng in text_lower:
            phase_state = state
            phase_confidence = "确认" if "confirmed" in state else "推定" if "presumed" in state else ""
            break

    return {
        "phenotypes": phenotypes,
        "zygosity": zygosity,
        "pathogenicity": pathogenicity,
        "phase_state": phase_state,
        "phase_confidence": phase_confidence,
    }


# ==================== 患者数量估算 ====================

def estimate_patient_count(text, sentences):
    """估算携带目标变异的患者数量。"""
    # 从匹配句子中查找数字
    patterns = [
        r'(\d+)\s*(?:patients?|cases?|probands?|subjects?|individuals?|carriers?)',
        r'(\d+)/\d+\s*(?:patients?|cases?)',
        r'n\s*=\s*(\d+)',
    ]
    for s in sentences:
        for p in patterns:
            m = re.search(p, s, re.IGNORECASE)
            if m:
                return int(m.group(1))
    return 1  # 默认 1


# ==================== 主搜索流程 ====================

def search_pdf(pdf_path, gene, keywords):
    """搜索单个 PDF 文件。返回结果字典。"""
    filename = os.path.basename(pdf_path)
    print(f"\n  搜索: {filename}")

    text, tables = extract_pdf_text(pdf_path)
    if not text:
        print(f"    -> 无法提取文本")
        return None

    print(f"    提取 {len(text)} 字符, {len(tables)} 个表格")

    # 搜索正文
    found, matched_kw, sentences = search_variant_in_text(text, keywords, gene)

    # 搜索表格
    matched_rows = search_variant_in_tables(tables, keywords, gene)

    # 合并表格中的发现
    for row_info in matched_rows:
        row_text = ' '.join(str(c) for c in row_info["row_content"] if c)
        if row_text not in sentences:
            sentences.append(row_text)
        if row_info["matched_keyword"] not in matched_kw:
            matched_kw.append(row_info["matched_keyword"])

    table_found = len(matched_rows) > 0
    total_found = found or table_found

    if total_found:
        print(f"    ✓ 提及目标变异 (正文: {found}, 表格: {table_found})")
        print(f"    匹配关键词: {matched_kw}")
    else:
        print(f"    - 未提及目标变异")

    # 提取表型信息
    phenotype_info = extract_phenotype_info(text, sentences if total_found else [])

    # 估算患者数
    patient_count = estimate_patient_count(text, sentences) if total_found else 0

    # 尝试提取标题 (前 200 字符中查找)
    title = text[:300].strip().replace('\n', ' ')[:200] if text else ""

    # 尝试获取 PMID
    pmid = ""
    pmid_match = re.search(r'PMID(\d+)', filename)
    if pmid_match:
        pmid = pmid_match.group(1)

    result = {
        "PMID": pmid,
        "文件": filename,
        "标题": title,
        "基因": gene,
        "cDNA变异": "",
        "蛋白变异": "",
        "变异提及": total_found,
        "匹配关键词": matched_kw,
        "正文匹配句": sentences[:10],
        "表格匹配行": [r["row_content"] for r in matched_rows],
        "患者数量": patient_count,
        "变异类型": "错义突变 (missense)" if total_found else "不适用",
        "致病性": ', '.join(phenotype_info["pathogenicity"]) if phenotype_info["pathogenicity"] else ("不适用" if not total_found else ""),
        "合子状态": ', '.join(phenotype_info["zygosity"]) if phenotype_info["zygosity"] else ("不适用" if not total_found else ""),
        "临床表型": ', '.join(phenotype_info["phenotypes"]) if phenotype_info["phenotypes"] else "",
        "相位状态": phenotype_info["phase_state"] if total_found else "",
        "相位置信度": phenotype_info["phase_confidence"],
        "表格数量": len(tables),
    }
    return result


# ==================== 主函数 ====================

def main():
    parser = argparse.ArgumentParser(
        description="离线 PDF 变异搜索器 — 仅在本地 PDF 中搜索目标变异，不联网",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
示例:
  # 基础用法
  python pdf_variant_search.py --pdf-dir D:\\project1\\sci --gene ABCA4 --variant "c.763C>T (p.Arg255Cys)"

  # 带转录本
  python pdf_variant_search.py --pdf-dir D:\\project1\\sci --gene ABCA4 --variant "c.1522T>C (p.Cys508Arg)" --transcript NM_001009944.3

  # 指定输出路径
  python pdf_variant_search.py --pdf-dir D:\\project1\\sci --gene ABCA4 --variant "c.763C>T (p.Arg255Cys)" --output results.json --excel-dir D:\\project1\\文献提取结果
        """,
    )
    parser.add_argument("--pdf-dir", required=True, help="PDF 文件所在目录路径")
    parser.add_argument("--gene", required=True, help="目标基因名称")
    parser.add_argument("--variant", required=True, help='目标变异，格式如 "c.763C>T (p.Arg255Cys)"')
    parser.add_argument("--transcript", default="", help="转录本 ID（可选，用于记录）")
    parser.add_argument("--output", default="pdf_variant_results.json", help="输出 JSON 路径")
    parser.add_argument("--excel-dir", default=r"D:\claude_code\project1\文献提取结果", help="Excel 汇总表输出目录")

    args = parser.parse_args()

    pdf_dir = args.pdf_dir
    if not os.path.isdir(pdf_dir):
        print(f"错误: PDF 目录不存在: {pdf_dir}")
        sys.exit(1)

    # 生成关键词
    keywords = generate_variant_keywords(args.variant)
    all_kw = []
    for cat in ["exact_cdna", "fuzzy_cdna", "protein"]:
        all_kw.extend(keywords.get(cat, []))

    print(f"目标变异: {args.variant}")
    print(f"基因: {args.gene}")
    if args.transcript:
        print(f"转录本: {args.transcript}")
    print(f"搜索关键词 ({len(all_kw)} 个):")
    print(f"  精确: {keywords['exact_cdna']}")
    print(f"  模糊: {keywords['fuzzy_cdna']}")
    print(f"  蛋白: {keywords['protein']}")
    print(f"\nPDF 目录: {pdf_dir}")

    # 查找所有 PDF 文件（包括 PMID*.pdf 和 PubMed*.pdf 和纯数字.pdf）
    pdf_files = []
    for f in sorted(os.listdir(pdf_dir)):
        if f.lower().endswith('.pdf'):
            pdf_files.append(os.path.join(pdf_dir, f))

    if not pdf_files:
        print(f"错误: 目录中无 PDF 文件: {pdf_dir}")
        sys.exit(1)

    print(f"找到 {len(pdf_files)} 个 PDF 文件\n")

    # 搜索所有 PDF
    results = []
    for pdf_path in pdf_files:
        result = search_pdf(pdf_path, args.gene, keywords)
        if result:
            results.append(result)

    # 统计
    mentioned = [r for r in results if r["变异提及"]]
    not_mentioned = [r for r in results if not r["变异提及"]]
    print(f"\n{'='*60}")
    print(f"搜索完成: {len(results)} 个 PDF")
    print(f"  提及目标变异: {len(mentioned)}")
    print(f"  未提及: {len(not_mentioned)}")

    for r in mentioned:
        print(f"\n  PMID {r['PMID']} ({r['文件']}):")
        print(f"    匹配关键词: {r['匹配关键词']}")
        print(f"    患者数: {r['患者数量']}")
        print(f"    致病性: {r['致病性']}")
        print(f"    合子状态: {r['合子状态']}")
        if r['正文匹配句']:
            print(f"    相关句: {r['正文匹配句'][0][:200]}...")

    # 保存 JSON
    with open(args.output, 'w', encoding='utf-8') as f:
        json.dump(results, f, ensure_ascii=False, indent=2)
    print(f"\nJSON 输出: {args.output}")

    # 生成 Excel CSV 汇总表
    os.makedirs(args.excel_dir, exist_ok=True)
    safe_variant = args.variant.replace('>', '_').replace(':', '_').replace('(', '').replace(')', '').replace(' ', '_')
    csv_path = os.path.join(args.excel_dir, f"{args.gene}_{safe_variant}_PDF文献汇总.csv")

    with open(csv_path, 'w', encoding='utf-8-sig', newline='') as f:
        writer = csv.writer(f)
        writer.writerow([
            "PMID", "标题", "是否提及此位点", "患者数", "致病性",
            "关联合子状态", "临床表型", "匹配关键词", "相关句子片段",
        ])
        for r in results:
            writer.writerow([
                r["PMID"],
                r["标题"][:200],
                "是" if r["变异提及"] else "否",
                r["患者数量"],
                r["致病性"],
                r["合子状态"],
                r["临床表型"],
                ', '.join(r["匹配关键词"]),
                (r["正文匹配句"][0] if r["正文匹配句"] else "")[:300],
            ])

    print(f"Excel 汇总表: {csv_path}")

    # 返回简要摘要
    print(f"\n{'='*60}")
    print("汇总:")
    print(f"{'PMID':<12} {'提及':<6} {'患者数':<6} {'致病性':<20} {'合子状态':<15}")
    print("-" * 60)
    for r in results:
        print(f"{r['PMID']:<12} {'是' if r['变异提及'] else '否':<6} {r['患者数量']:<6} {r['致病性'][:18]:<20} {r['合子状态'][:13]:<15}")


if __name__ == "__main__":
    main()