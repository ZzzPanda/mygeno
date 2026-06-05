#!/usr/bin/env python3
"""
PubMed Parser MCP Server
使用 FastMCP 框架提供 PubMed 文献查询和变异信息提取服务。

Usage:
  python server.py
  # 或直接运行: fastmcp run server:server
"""

from fastmcp import FastMCP

# 导入 pubmed_client 核心功能
from pubmed_client import (
    fetch_ncbi_abstract,
    fetch_europe_pmc_fulltext,
    fetch_europe_pmc_text,
    fetch_pmc_fulltext,
    parse_ncbi_xml,
    parse_pmc_europe_xml,
    build_variant_keywords,
    find_variant_sentences,
    split_sentences,
    extract_pathogenicity,
    extract_zygosity,
    extract_inheritance,
    extract_patient_phenotypes,
    extract_clinical_details,
    extract_co_variants,
    infer_variant_type,
    expand_protein_keywords,
)

# 创建 FastMCP 服务器实例
mcp = FastMCP("PubMed Parser")

# ── 工具函数 ──────────────────────────────────────────────────────────────────


@mcp.tool()
def search_pubmed_by_pmid(pmids: list[str]) -> list[dict]:
    """
    根据 PMID 列表获取 PubMed 文章摘要信息。

    Args:
        pmids: PubMed ID 列表，如 ["12345678", "23456789"]

    Returns:
        文章信息列表，包含标题、摘要、作者、期刊、年份等
    """
    results = []
    for pmid in pmids:
        pmid = pmid.strip()
        if not pmid:
            continue

        # 尝试多种数据源
        xml_text = fetch_ncbi_abstract(pmid)
        if xml_text:
            result = parse_ncbi_xml(xml_text)
            if result and result.get("PMID"):
                results.append(result)
                continue

        # Europe PMC fallback
        europe_data = fetch_europe_pmc_text(pmid)
        if europe_data:
            results.append({
                "PMID": pmid,
                "标题": europe_data.get("title", ""),
                "摘要": europe_data.get("abstract", ""),
                "作者": europe_data.get("authors", []),
                "期刊": europe_data.get("journal", ""),
                "发表年份": europe_data.get("year", ""),
                "全文来源": "europe_pmc",
            })
            continue

        results.append({"PMID": pmid, "错误": "未找到文章"})

    return results


@mcp.tool()
def get_fulltext_by_pmid(pmid: str) -> dict:
    """
    获取 PubMed 文章的全文内容（优先 PMC 开放获取全文）。

    Args:
        pmid: PubMed ID

    Returns:
        包含全文内容的字典，包含 '全文'、'全文来源' 等字段
    """
    # 优先尝试 PMC 全文
    xml_text = fetch_pmc_fulltext(pmid)
    if xml_text:
        result = parse_pmc_europe_xml(xml_text, {"PMID": pmid})
        if result.get("全文"):
            return result

    # Europe PMC 全文 fallback
    xml_text = fetch_europe_pmc_fulltext(pmid)
    if xml_text:
        result = parse_pmc_europe_xml(xml_text, {"PMID": pmid})
        if result.get("全文"):
            return result

    # 最终尝试 NCBI 摘要
    xml_text = fetch_ncbi_abstract(pmid)
    if xml_text:
        result = parse_ncbi_xml(xml_text)
        if result and result.get("PMID"):
            return result

    return {"PMID": pmid, "错误": "无法获取全文", "全文": "", "全文来源": "none"}


@mcp.tool()
def extract_variant_info(
    pmid: str,
    cdna: str = "",
    protein: str = "",
    gene: str = "",
    transcript: str = "",
) -> dict:
    """
    从 PubMed 文章中提取目标变异的信息。

    Args:
        pmid: PubMed ID
        cdna: cDNA 变异，如 "c.1166G>A"
        protein: 蛋白变异，如 "p.R389H" 或 "R389H"
        gene: 基因符号，如 "ABCG5"
        transcript: 转录本 ID，如 "NM_022436.3"

    Returns:
        变异信息字典，包含匹配句子、致病性、遗传方式、患者表型等
    """
    # 获取全文
    xml_text = fetch_pmc_fulltext(pmid)
    result = {}
    if xml_text:
        result = parse_pmc_europe_xml(xml_text, {"PMID": pmid})

    if not result.get("全文"):
        xml_text = fetch_europe_pmc_fulltext(pmid)
        if xml_text:
            result = parse_pmc_europe_xml(xml_text, {"PMID": pmid})

    if not result.get("全文"):
        xml_text = fetch_ncbi_abstract(pmid)
        if xml_text:
            result = parse_ncbi_xml(xml_text)

    if not result or not result.get("全文"):
        return {"PMID": pmid, "错误": "无法获取文章内容"}

    # 构建关键词
    keywords = build_variant_keywords(cdna, protein, transcript)

    # 查找变异提及句子
    variant_sentences, matched_kws = find_variant_sentences(result.get("全文", ""), keywords)

    # 提取各项信息
    full_text_lower = result.get("全文", "").lower()

    info = {
        "PMID": pmid,
        "标题": result.get("标题", ""),
        "变异提及": bool(variant_sentences),
        "匹配关键词": matched_kws,
        "相关句子": variant_sentences[:20],  # 最多返回 20 条
        "变异类型": infer_variant_type(cdna, protein, variant_sentences),
        "致病性": extract_pathogenicity(variant_sentences, full_text_lower),
        "合子状态": extract_zygosity(
            variant_sentences,
            full_text_lower,
            target_cdna=cdna,
            target_protein=protein,
            keywords=keywords,
        ),
        "遗传模式": extract_inheritance(variant_sentences, full_text_lower),
        "临床表型": extract_patient_phenotypes(variant_sentences, keywords=keywords),
        "临床详情": extract_clinical_details(variant_sentences, keywords=keywords),
    }

    # 提取共存变异
    co_variants = extract_co_variants(
        variant_sentences,
        target_cdna=cdna,
        target_protein=protein,
        variant_keywords=keywords,
    )
    if co_variants:
        info["共存变异"] = co_variants

    return info


@mcp.tool()
def search_variant_keywords(protein: str = "", cdna: str = "") -> dict:
    """
    生成变异搜索关键词的所有变体（用于调试或在其他工具中使用）。

    Args:
        protein: 蛋白变异，如 "p.R389H" 或 "R389H"
        cdna: cDNA 变异，如 "c.1166G>A"

    Returns:
        关键词变体字典
    """
    keywords = build_variant_keywords(cdna, protein)
    return {
        "原始输入": {"protein": protein, "cdna": cdna},
        "精确匹配关键词": keywords["exact"],
        "模糊匹配关键词": keywords["fuzzy"],
        "蛋白变体关键词": list(keywords["protein"]),
        "历史命名": keywords["historical"],
    }


@mcp.tool()
def batch_extract_variants(
    pmids: list[str],
    cdna: str = "",
    protein: str = "",
    gene: str = "",
    transcript: str = "",
) -> dict:
    """
    批量从多个 PubMed 文章中提取变异信息。

    Args:
        pmids: PubMed ID 列表
        cdna: cDNA 变异
        protein: 蛋白变异
        gene: 基因符号
        transcript: 转录本 ID

    Returns:
        包含所有文章结果的字典
    """
    results = []
    for pmid in pmids:
        pmid = pmid.strip()
        if not pmid:
            continue
        try:
            info = extract_variant_info(pmid, cdna, protein, gene, transcript)
            results.append(info)
        except Exception as e:
            results.append({"PMID": pmid, "错误": str(e)})

    return {
        "查询参数": {"pmids": pmids, "cdna": cdna, "protein": protein, "gene": gene},
        "结果数量": len(results),
        "结果": results,
    }


# ── 资源 ──────────────────────────────────────────────────────────────────────


@mcp.resource("pubmed://pmid/{pmid}")
def pubmed_article(pmid: str) -> str:
    """获取指定 PMID 的 PubMed 文章信息（作为文本资源）。"""
    articles = search_pubmed_by_pmid([pmid])
    if articles:
        a = articles[0]
        lines = [
            f"PMID: {a.get('PMID', '')}",
            f"标题: {a.get('标题', '')}",
            f"期刊: {a.get('期刊', '')} ({a.get('发表年份', '')})",
            f"作者: {', '.join(a.get('作者', [])[:5])}",
            "",
            "摘要:",
            a.get('摘要', '无摘要'),
        ]
        return "\n".join(lines)
    return f"未找到 PMID: {pmid}"


# ── 提示 (Prompts) ────────────────────────────────────────────────────────────


@mcp.prompt()
def analyze_patient_variants(pmid: str, gene: str) -> str:
    """生成分析某基因相关患者变异的提示词。"""
    return f"""分析以下 PubMed 文章中 {gene} 基因相关患者的变异信息。

PMID: {pmid}

请提取:
1. 患者的变异类型（cDNA 和蛋白改变）
2. 致病性评级
3. 遗传模式（显性/隐性）
4. 患者表型（临床特征）
5. 是否为复合杂合或纯合变异
"""


if __name__ == "__main__":
    mcp.run()