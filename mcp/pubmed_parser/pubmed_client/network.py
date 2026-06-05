"""
Network fetching functions for PubMed/PMC/Europe PMC APIs.
"""

import json
import time
import random
import urllib.request
import urllib.parse
import urllib.error

from .constants import (
    MAX_RETRIES,
    DAILY_SITE_LIMIT,
    check_and_increment,
)


# ── 网络请求 ──

def safe_request(url, description="", retries=MAX_RETRIES):
    req = urllib.request.Request(url)
    req.add_header("User-Agent", "PubMedVariantExtractor/7.0")
    if description:
        print(f"  -> {description}")
    for attempt in range(retries + 1):
        try:
            if attempt > 0:
                wait = random.randint(20, 60) * attempt
                print(f"  重试 {attempt}/{retries}，等待 {wait}s...")
                time.sleep(wait)
            with urllib.request.urlopen(req, timeout=45) as response:
                return response.read().decode("utf-8", errors="replace")
        except urllib.error.HTTPError as e:
            if e.code == 429:
                time.sleep(random.randint(30, 90))
                continue
        except Exception:
            pass
        if attempt < retries:
            time.sleep(random.randint(10, 40))
    return None


def fetch_ncbi_abstract(pmid):
    site = "eutils.ncbi.nlm.nih.gov"
    allowed, count = check_and_increment(site)
    if not allowed:
        print(f"\n  [警告] {site} 今日已访问 {count} 次（上限 {DAILY_SITE_LIMIT}）")
        return None
    params = {"db": "pubmed", "id": pmid, "rettype": "xml", "retmode": "text"}
    url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?{urllib.parse.urlencode(params)}"
    return safe_request(url, f"NCBI abstract PMID:{pmid}")


def fetch_europe_pmc_fulltext(pmid):
    site = "www.ebi.ac.uk"
    allowed, count = check_and_increment(site)
    if not allowed:
        print(f"\n  [警告] {site} 今日已访问 {count} 次")
        return None
    check_url = (f"https://www.ebi.ac.uk/europepmc/webservices/rest/"
                 f"search?query=EXT_ID:{pmid}&format=json&resultType=core")
    data = safe_request(check_url, f"Europe PMC check PMID:{pmid}")
    if not data:
        return None
    try:
        result = json.loads(data)
        entries = result.get("resultList", {}).get("result", [])
        if not entries:
            return None
        entry = entries[0]
        if entry.get("isOpenAccess") == "Y" and entry.get("pmcid"):
            pmcid = entry["pmcid"]
            allowed2, _ = check_and_increment(site)
            if not allowed2:
                print(f"  [警告] {site} 今日已达上限")
                return None
            fulltext_url = (f"https://www.ebi.ac.uk/europepmc/webservices/rest/"
                           f"{pmcid}/fullTextXML")
            return safe_request(fulltext_url, f"Europe PMC XML {pmcid}")
    except (json.JSONDecodeError, KeyError):
        pass
    return None


def fetch_pmc_fulltext(pmid):
    """通过 PMC esearch 获取 PMC 全文 XML（通过 PMID 搜索对应 PMC 文章）"""
    site = "eutils.ncbi.nlm.nih.gov"
    allowed, count = check_and_increment(site)
    if not allowed:
        print(f"\n  [警告] {site} 今日已访问 {count} 次")
        return None
    # 使用 esearch 在 PMC 中搜索对应 PMID 的文章（比 elink 更准确）
    search_url = (f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
                  f"?db=pmc&term={pmid}[pmid]&retmode=json")
    data = safe_request(search_url, f"PMC esearch PMID:{pmid}")
    if not data:
        return None
    try:
        result = json.loads(data)
        pmc_ids = result.get("esearchresult", {}).get("idlist", [])
        if not pmc_ids:
            return None
        xml_url = (f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
                   f"?db=pmc&id={pmc_ids[0]}&rettype=xml&retmode=text")
        return safe_request(xml_url, f"PMC XML PMC{pmc_ids[0]}")
    except (json.JSONDecodeError, KeyError, IndexError):
        pass
    return None


# ── 转录本编号转换查询 ──

def fetch_ncbi_gene_info(transcript_id):
    """
    通过 NCBI Gene/RefSeq 接口查询转录本信息，
    返回该转录本对应的 CDS 起始位置等元数据，
    用于计算不同转录本版本间的 cDNA 编号偏移量。
    如果查询失败返回 None。
    """
    if not transcript_id:
        return None
    site = "eutils.ncbi.nlm.nih.gov"
    allowed, _ = check_and_increment(site)
    if not allowed:
        return None

    # 用 efetch 获取转录本的 FASTA + 元数据
    url = (f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
           f"?db=nucleotide&id={transcript_id}&rettype=fasta&retmode=text")
    data = safe_request(url, f"NCBI nucleotide {transcript_id}", retries=1)
    if not data:
        return None

    # 从 FASTA header 中提取 CDS 信息
    # Header 格式类似: >ref|NM_022436.3| Homo sapiens ATP binding cassette subfamily G member 5 (ABCG5), ...
    header_line = data.split('\n')[0] if data else ""
    info = {"transcript": transcript_id, "header": header_line}

    # 尝试用 esummary 获取更多结构化信息
    summary_url = (f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
                   f"?db=nucleotide&id={transcript_id}&retmode=json")
    summary_data = safe_request(summary_url, f"NCBI summary {transcript_id}", retries=1)
    if summary_data:
        try:
            sresult = json.loads(summary_data)
            uid = list(sresult.get("result", {}).keys())
            if uid and uid[0] != "uids":
                summary_info = sresult["result"][uid[0]]
                info["title"] = summary_info.get("title", "")
                info["extra"] = summary_info.get("extra", "")
        except (json.JSONDecodeError, KeyError):
            pass

    return info


def fetch_europe_pmc_text(pmid):
    """
    通过 Europe PMC 搜索 API 获取文章的摘要文本（作为 fallback）。
    当 NCBI 摘要和 PMC 全文都不可用时使用。
    返回 dict 包含标题、摘要、作者等信息；失败返回 None。
    """
    site = "www.ebi.ac.uk"
    allowed, count = check_and_increment(site)
    if not allowed:
        print(f"  [警告] {site} 今日已访问 {count} 次")
        return None
    url = (f"https://www.ebi.ac.uk/europepmc/webservices/rest/"
           f"search?query=EXT_ID:{pmid}&format=json&resultType=core")
    data = safe_request(url, f"Europe PMC text PMID:{pmid}", retries=1)
    if not data:
        return None
    try:
        result = json.loads(data)
        entries = result.get("resultList", {}).get("result", [])
        if not entries:
            return None
        entry = entries[0]
        pub_info = {
            "title": entry.get("title", ""),
            "abstract": entry.get("abstractText", ""),
            "journal": entry.get("journalTitle", ""),
            "year": str(entry.get("pubYear", "")),
            "authors": [],
            "pmid": pmid,
            "pmcid": entry.get("pmcid", ""),
            "isOpenAccess": entry.get("isOpenAccess", "N"),
        }
        author_str = entry.get("authorString", "")
        if author_str:
            pub_info["authors"] = [a.strip() for a in author_str.split(",") if a.strip()]
        return pub_info
    except (json.JSONDecodeError, KeyError):
        pass
    return None


def fetch_pmc_html(pmid):
    """
    从 PMC 网页获取文章的 HTML 全文（用于出版商限制 XML 时的 fallback）。
    PMC 网站对大部分文章提供免费 HTML 全文，即使 XML 下载受限。
    返回 HTML 文本；失败返回 None。
    """
    site = "www.ncbi.nlm.nih.gov"
    allowed, count = check_and_increment(site)
    if not allowed:
        print(f"  [警告] {site} 今日已访问 {count} 次")
        return None
    # 先通过 esearch 获取 PMCID
    search_url = (f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
                  f"?db=pmc&term={pmid}[pmid]&retmode=json")
    data = safe_request(search_url, f"PMC esearch HTML PMID:{pmid}", retries=1)
    if not data:
        return None
    try:
        result = json.loads(data)
        pmc_ids = result.get("esearchresult", {}).get("idlist", [])
        if not pmc_ids:
            return None
        pmcid = pmc_ids[0]
    except (json.JSONDecodeError, KeyError, IndexError):
        return None

    # 获取 PMC HTML 页面
    html_url = f"https://www.ncbi.nlm.nih.gov/pmc/articles/{pmcid}/"
    allowed2, _ = check_and_increment(site)
    if not allowed2:
        return None
    html_data = safe_request(html_url, f"PMC HTML {pmcid}", retries=1)
    return html_data