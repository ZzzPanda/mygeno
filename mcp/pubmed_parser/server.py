#!/usr/bin/env python3
"""
PubMed Parser MCP Server
使用 FastMCP 框架提供 PDF 文献解析和变异信息提取服务。

Usage:
  python server.py
  # 或直接运行: fastmcp run server:server
"""

from fastmcp import FastMCP

import pymupdf

# 导入 pubmed_client 核心功能
from pubmed_client import (
    extract_text_from_pdf,
    extract_tables_from_pdf,
    build_variant_keywords,
    find_variant_sentences,
    split_sentences,
    extract_pathogenicity,
    extract_zygosity,
    extract_inheritance,
    extract_patient_phenotypes,
    extract_clinical_details,
    extract_co_variants,
    extract_phase_evidence,
    extract_patient_count,
    infer_variant_type,
)

# 创建 FastMCP 服务器实例
mcp = FastMCP("PubMed Parser")


def _build_result_dict(pdf_path, title, gene, cdna, protein, transcript,
                       full_text, tables, variant_sentences, matched_kws,
                       table_rows, table_count):
    """构建符合 target 格式的结果字典。"""
    full_text_lower = full_text.lower()
    keywords = build_variant_keywords(cdna, protein, transcript)

    # 变异类型
    variant_type = infer_variant_type(cdna, protein, variant_sentences)

    # 致病性
    pathogenicity = extract_pathogenicity(variant_sentences, full_text_lower)

    # 合子状态（传入tables和keywords）
    zygosity = extract_zygosity(
        variant_sentences,
        full_text_lower,
        target_cdna=cdna,
        target_protein=protein,
        tables=tables,
        keywords=keywords,
    )

    # 遗传模式
    inheritance = extract_inheritance(variant_sentences, full_text_lower)

    # 临床表型
    phenotypes = extract_patient_phenotypes(variant_sentences, tables, keywords)
    phenotype_str = "、".join(phenotypes) if phenotypes else ""

    # 临床详情
    clinical_details = extract_clinical_details(variant_sentences, tables, keywords)

    # 相位信息
    co_variants = extract_co_variants(
        variant_sentences,
        tables=tables,
        target_cdna=cdna,
        target_protein=protein,
        variant_keywords=keywords,
    )
    phase_result, _ = extract_phase_evidence(variant_sentences, co_variants, zygosity)

    # 患者数量
    patient_count = extract_patient_count(variant_sentences, tables, keywords, cdna)

    # 提取PMID从文件名
    pmid = ""
    if pdf_path:
        import re
        m = re.search(r'PMID:(\d+)', pdf_path)
        if m:
            pmid = m.group(1)

    return {
        "PMID": pmid,
        "文件": pdf_path.split("/")[-1] if pdf_path else "",
        "标题": title,
        "基因": gene.upper() if gene else "",
        "cDNA变异": cdna,
        "蛋白变异": protein,
        "变异提及": bool(variant_sentences),
        "匹配关键词": matched_kws,
        "正文匹配句": variant_sentences[:20],
        "表格匹配行": table_rows[:20],
        "患者数量": patient_count,
        "变异类型": variant_type,
        "致病性": pathogenicity,
        "合子状态": zygosity,
        "临床表型": phenotype_str,
        "相位状态": phase_result.get("phase_status", ""),
        "相位置信度": phase_result.get("phase_confidence", ""),
        "表格数量": table_count,
    }


@mcp.tool()
def extract_variant_info(
    pdf_path: str,
    cdna: str = "",
    protein: str = "",
    gene: str = "",
    transcript: str = "",
) -> dict:
    """
    从 PDF 文件中提取目标变异的信息。

    Args:
        pdf_path: PDF 文件路径（如 "/Users/roger/Documents/GitHub/mygeno/data/pdf/PMID:33374015.pdf"）
        cdna: cDNA 变异，如 "c.2279C>T"
        protein: 蛋白变异，如 "p.Thr760Met" 或 "Thr760Met"
        gene: 基因符号，如 "CFTR"
        transcript: 转录本 ID，如 "NM_000492.4"

    Returns:
        变异信息字典，包含匹配句子、致病性、遗传方式、患者表型等
    """
    # 1. 提取文本
    try:
        full_text = extract_text_from_pdf(pdf_path)
    except Exception as e:
        return {"错误": f"无法读取 PDF: {e}"}

    if not full_text:
        return {"错误": "PDF 内容为空"}

    # 2. 提取表格
    tables = extract_tables_from_pdf(pdf_path)

    # 3. 构建关键词
    keywords = build_variant_keywords(cdna, protein, transcript)

    # 4. 查找变异提及句子
    variant_sentences, matched_kws = find_variant_sentences(full_text, keywords)

    # 5. 收集表格中包含变异的行
    table_rows = []
    for table_info in tables:
        for row in table_info["rows"]:
            row_text = "\t".join(str(cell or "") for cell in row)
            row_lower = row_text.lower()
            for kw in keywords.get("all", []):
                if kw.lower() in row_lower:
                    table_rows.append(row)
                    break

    # 6. 提取标题（从 PDF 元数据或第一行）
    title = ""
    try:
        import pymupdf
        with pymupdf.open(pdf_path) as doc:
            metadata = doc.metadata
            if metadata.get("title"):
                title = metadata["title"]
    except Exception:
        pass
    if not title and variant_sentences:
        # 尝试从第一句话提取
        title = variant_sentences[0][:200] if variant_sentences else ""

    # 7. 构建结果
    return _build_result_dict(
        pdf_path=pdf_path,
        title=title,
        gene=gene,
        cdna=cdna,
        protein=protein,
        transcript=transcript,
        full_text=full_text,
        tables=tables,
        variant_sentences=variant_sentences,
        matched_kws=matched_kws,
        table_rows=table_rows,
        table_count=len(tables),
    )


@mcp.tool()
def parse_pdf(pdf_path: str) -> dict:
    """
    解析 PDF 文件，提取文本和表格。

    Args:
        pdf_path: PDF 文件路径

    Returns:
        包含 'text', 'tables', 'table_count' 的字典
    """
    try:
        text = extract_text_from_pdf(pdf_path)
    except Exception as e:
        return {"错误": f"无法读取 PDF: {e}", "text": "", "tables": [], "table_count": 0}

    tables = extract_tables_from_pdf(pdf_path)

    return {
        "text": text,
        "tables": tables,
        "table_count": len(tables),
    }


@mcp.tool()
def analyze_variant(
    text: str,
    cdna: str = "",
    protein: str = "",
    gene: str = "",
    transcript: str = "",
) -> dict:
    """
    从文本（论文文本/摘要）中提取变异信息。

    Args:
        text: 论文文本或摘要
        cdna: cDNA 变异，如 "c.2279C>T"
        protein: 蛋白变异，如 "p.Thr760Met" 或 "Thr760Met"
        gene: 基因符号，如 "CFTR"
        transcript: 转录本 ID，如 "NM_000492.4"

    Returns:
        变异信息字典
    """
    if not text:
        return {"错误": "文本为空"}

    # 构建关键词
    keywords = build_variant_keywords(cdna, protein, transcript)

    # 查找变异提及句子
    variant_sentences, matched_kws = find_variant_sentences(text, keywords)

    # 表格为空（文本输入不解析表格）
    tables = []

    # 收集表格匹配行（文本输入无表格）
    table_rows = []

    # 标题从文本提取
    title = text[:200] if text else ""

    return _build_result_dict(
        pdf_path="",
        title=title,
        gene=gene,
        cdna=cdna,
        protein=protein,
        transcript=transcript,
        full_text=text,
        tables=tables,
        variant_sentences=variant_sentences,
        matched_kws=matched_kws,
        table_rows=table_rows,
        table_count=0,
    )


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


# ── 资源 ──────────────────────────────────────────────────────────────────────


@mcp.resource("pdf://parse/{pdf_path}")
def parse_pdf_resource(pdf_path: str) -> str:
    """解析指定路径的 PDF 文件并返回文本内容。"""
    try:
        result = parse_pdf(pdf_path)
        if "错误" in result:
            return f"错误: {result['错误']}"
        lines = [
            f"文件: {pdf_path}",
            f"表格数量: {result['table_count']}",
            "",
            "文本内容（前 2000 字符）:",
            result["text"][:2000],
        ]
        return "\n".join(lines)
    except Exception as e:
        return f"解析失败: {e}"


# ── 提示 (Prompts) ────────────────────────────────────────────────────────────


@mcp.prompt()
def analyze_patient_variants(pdf_path: str, gene: str) -> str:
    """生成分析某基因相关患者变异的提示词。"""
    return f"""分析以下 PDF 文献中 {gene} 基因相关患者的变异信息。

文件路径: {pdf_path}

请提取:
1. 患者的变异类型（cDNA 和蛋白改变）
2. 致病性评级
3. 遗传模式（显性/隐性）
4. 患者表型（临床特征）
5. 是否为复合杂合或纯合变异
6. 相位状态（顺式/反式）
"""


if __name__ == "__main__":
    import sys
    if "--test" in sys.argv:
        # 测试模式：不启动 MCP 服务器，只验证导入和简单调用
        print("Running test mode...")
        from pubmed_client import (
            extract_text_from_pdf,
            extract_tables_from_pdf,
            build_variant_keywords,
            find_variant_sentences,
        )
        print("pubmed_client import successful")

        # 测试关键词构建
        keywords = build_variant_keywords("c.2279C>T", "p.Thr760Met")
        print(f"Keywords built: {list(keywords['all'])[:5]}")

        print("\nTest complete, server not started.")
        print("To start MCP server, connect via Claude Code MCP configuration.")
        sys.exit(0)

    mcp.run()