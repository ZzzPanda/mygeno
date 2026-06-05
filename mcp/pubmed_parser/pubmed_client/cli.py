"""
CLI interface for PubMed variant extraction.
"""

import argparse
import json
import time
import re
import os

from .constants import DAILY_SITE_LIMIT, load_site_counts
from .network import (
    fetch_ncbi_abstract,
    fetch_europe_pmc_fulltext,
    fetch_europe_pmc_text,
    fetch_pmc_fulltext,
    fetch_pmc_html,
    fetch_ncbi_gene_info,
    safe_request,
)
from .constants import random_sleep, periodic_long_pause
from .xml_parser import parse_ncbi_xml, parse_pmc_europe_xml, parse_pmc_html, _empty_result
from .variants import build_variant_keywords
from .extractors import extract_info_for_variant
from .summarizer import (
    generate_summary_paragraph,
    generate_one_sentence_summary,
    _generate_excel_csv,
)


def main():
    parser = argparse.ArgumentParser(
        description="从PubMed/PMC文献中仅针对目标变异提取信息"
    )
    parser.add_argument("--pmids", nargs="+", required=True, help="PubMed ID列表")
    parser.add_argument("--gene", default=None, help="目标基因名称")
    parser.add_argument("--variant", default=None,
                        help='目标变异，如 "c.1166G>A (p.Arg389His)"')
    parser.add_argument("--transcript", default=None,
                        help="转录本ID，如 NM_022436.3")
    parser.add_argument("--output", default="pubmed_variant_results.json",
                        help="输出JSON文件路径")
    parser.add_argument("--excel-dir", default=r"./pubmed_extractor_output",
                        help="Excel汇总表输出目录，默认 ./pubmed_extractor_output")
    parser.add_argument("--pdf-results", default=None,
                        help="PDF搜索结果的JSON文件路径（用于交叉引用，补充在线API遗漏的变异检出）")
    args = parser.parse_args()

    # 解析目标变异
    target_cdna = ""
    target_protein = ""
    if args.variant:
        cdna_m = re.search(r'c\.([\d\w_*><+=\-\+]+)', args.variant)
        if cdna_m:
            target_cdna = f"c.{cdna_m.group(1)}"
        prot_m = re.search(r'p\.\s*([\w?*]+)', args.variant)
        if prot_m:
            target_protein = f"p.{prot_m.group(1)}"

    if not args.gene:
        print("提示: 建议提供 --gene 参数以精确定位目标基因")

    # 去重
    seen = set()
    unique_pmids = []
    for pmid in args.pmids:
        if pmid not in seen:
            seen.add(pmid)
            unique_pmids.append(pmid)
    if len(unique_pmids) < len(args.pmids):
        print(f"去重: {len(args.pmids)} -> {len(unique_pmids)} 篇")

    # 构建变异关键词
    keywords = build_variant_keywords(target_cdna, target_protein, args.transcript or "")
    print(f"\n目标变异: {target_cdna} {target_protein}")
    print(f"搜索关键词 ({len(keywords['all'])} 个):")
    print(f"  精确: {keywords['exact']}")
    print(f"  模糊: {keywords['fuzzy']}")
    print(f"  蛋白: {keywords['protein']}")

    # 转录本信息查询
    if args.transcript:
        print(f"\n查询转录本信息: {args.transcript}")
        tx_info = fetch_ncbi_gene_info(args.transcript)
        if tx_info:
            print(f"  转录本信息: {tx_info.get('title', '')}")
        time.sleep(5)

    # 站点访问计数状态
    counts = load_site_counts()
    print(f"\n今日 ({__import__('datetime').date.today()}) 站点访问计数:")
    for site_key, count in counts.items():
        status = "[已达上限]" if count >= DAILY_SITE_LIMIT else "[OK]"
        print(f"   {site_key}: {count}/{DAILY_SITE_LIMIT} {status}")

    # v9: 加载PDF搜索结果用于交叉引用
    pdf_data = {}
    if args.pdf_results and os.path.exists(args.pdf_results):
        with open(args.pdf_results, 'r', encoding='utf-8') as f:
            pdf_list = json.load(f)
            for pr in pdf_list:
                pid = str(pr.get("PMID", "")).strip()
                # 回退: 从文件名提取 PMID (如 PubMed22326530.pdf → 22326530)
                if not pid:
                    fname = pr.get("文件", "")
                    m = re.search(r'(\d{7,9})', fname)
                    if m:
                        pid = m.group(1)
                if pid:
                    # 如果同一个 PMID 有多个条目（文件名不同），保留提及变异的
                    if pid in pdf_data and not pr.get("变异提及"):
                        continue
                    pdf_data[pid] = pr
        print(f"\nPDF交叉引用: 已加载 {len(pdf_data)} 条PDF搜索结果")
    elif args.pdf_results:
        print(f"\n[警告] PDF结果文件不存在: {args.pdf_results}")

    results = []
    request_count = 0
    stopped_early = False

    for i, pmid in enumerate(unique_pmids):
        print(f"\n{'='*80}")
        print(f"PMID: {pmid}  ({i+1}/{len(unique_pmids)})")
        print(f"{'='*80}")

        # Step 1: 获取NCBI摘要
        xml_text = fetch_ncbi_abstract(pmid)
        if xml_text is None and load_site_counts().get("eutils.ncbi.nlm.nih.gov", 0) >= DAILY_SITE_LIMIT:
            stopped_early = True
            break
        request_count += 1
        periodic_long_pause(request_count)

        result = None
        if xml_text:
            result = parse_ncbi_xml(xml_text)

        # NCBI 摘要获取失败时，回退到 Europe PMC
        if not result:
            print(f"  NCBI 摘要获取失败，尝试 Europe PMC...")
            epmc_info = fetch_europe_pmc_text(pmid)
            if epmc_info:
                result = _empty_result()
                result["PMID"] = pmid
                result["标题"] = epmc_info.get("title", "")
                result["摘要"] = epmc_info.get("abstract", "")
                result["全文"] = epmc_info.get("abstract", "")
                result["全文来源"] = "epmc_abstract"
                result["作者"] = epmc_info.get("authors", [])
                result["期刊"] = epmc_info.get("journal", "")
                result["发表年份"] = epmc_info.get("year", "")
                print(f"  已从 Europe PMC 获取摘要")
            else:
                print(f"  无法获取任何来源，跳过")
                continue

        # Step 2: 随机间隔
        random_sleep()
        request_count += 1
        periodic_long_pause(request_count)

        # Step 3: 尝试获取全文 (Europe PMC -> PMC)
        fulltext_xml = fetch_europe_pmc_fulltext(pmid)
        request_count += 1
        if not fulltext_xml:
            random_sleep()
            request_count += 1
            fulltext_xml = fetch_pmc_fulltext(pmid)

        if fulltext_xml:
            result = parse_pmc_europe_xml(fulltext_xml, result)
            print(f"  全文获取成功 (来源: {result['全文来源']})")
            if result.get("tables"):
                print(f"  提取到 {len(result['tables'])} 个表格")
            # 如果 PMC XML 受限，尝试从 PMC HTML 页面获取全文（网页版通常免费）
            if result.get("全文来源") == "pmc_restricted":
                print(f"  PMC XML 受出版商限制，尝试 PMC HTML 网页版...")
                html_data = fetch_pmc_html(pmid)
                if html_data:
                    result = parse_pmc_html(html_data, result)
                    print(f"  已获取 PMC HTML (来源: {result['全文来源']})")
                    if result.get("tables"):
                        print(f"  从 HTML 提取到 {len(result['tables'])} 个表格")
                else:
                    # 回退到 Europe PMC 摘要
                    print(f"  PMC HTML 获取失败，尝试 Europe PMC 摘要...")
                    epmc_text = fetch_europe_pmc_text(pmid)
                    if epmc_text and epmc_text.get("abstract"):
                        result["全文"] = epmc_text.get("abstract", "")
                        result["全文来源"] = "epmc_abstract"
                        print(f"  已获取 Europe PMC 摘要文本 ({len(result['全文'])} 字符)")
        else:
            print(f"  仅获取摘要")
            # NCBI 摘要成功但无全文，尝试补充 Europe PMC 摘要
            epmc_text = fetch_europe_pmc_text(pmid)
            if epmc_text and epmc_text.get("abstract"):
                result["全文"] = epmc_text.get("abstract", "")
                result["全文来源"] = "epmc_abstract"
                print(f"  已补充 Europe PMC 摘要文本 ({len(result['全文'])} 字符)")

        # Step 4: 仅针对目标变异提取信息
        pdf_match = pdf_data.get(str(pmid))
        result = extract_info_for_variant(result, args.gene, target_cdna, target_protein, keywords, pdf_match)

        # Step 5: 生成一句话概括和总结段落
        result["一句话概括"] = generate_one_sentence_summary(result)
        if result["变异提及"]:
            result["总结段落"] = generate_summary_paragraph(result)
        else:
            result["总结段落"] = result["一句话概括"]  # 未提及时使用文献简介作为总结

        # 打印结果
        print(f"\n  标题: {result['标题']}")
        print(f"  一句话概括: {result['一句话概括']}")
        print(f"  作者: {', '.join(result['作者'][:3])}{'...' if len(result['作者']) > 3 else ''}")
        print(f"  期刊: {result['期刊']} ({result['发表年份']})")
        print(f"  全文来源: {result['全文来源']}")

        if result["变异提及"]:
            print(f"\n  目标变异: {result['基因']} {result['cDNA变异']} {result['蛋白变异']}")
            print(f"    变异类型: {result['变异类型']}")
            print(f"    致病性:   {result['致病性']}")
            print(f"    遗传方式: {result['合子状态']}")
            print(f"    临床表型: {result['临床表型']}")
            print(f"    遗传模式: {result['遗传模式']}")
            print(f"    匹配关键词: {', '.join(result.get('匹配关键词', []))}")

            if result.get("共存变异"):
                co_var_display = []
                for cv in result['共存变异']:
                    if isinstance(cv, dict):
                        parts = [cv.get("cdna") or "", cv.get("蛋白变异") or ""]
                        co_var_display.append(" ".join(p for p in parts if p))
                    else:
                        co_var_display.append(str(cv))
                if co_var_display:
                    print(f"    共存变异: {', '.join(co_var_display)}")
            if result.get("反式确认"):
                print(f"    反式位置确认: 是")
            if result.get("顺式确认"):
                print(f"    顺式位置确认: 是")
            if result.get("相位状态"):
                print(f"    相位状态: {result['相位状态']} ({result.get('相位置信度', '')})")
                print(f"    相位详情: {result.get('相位详情', '')}")
            if result.get("亲本检测"):
                print(f"    亲本检测: 已进行")
                if result.get("母源变异"):
                    print(f"      母源变异: {result['母源变异']}")
                if result.get("父源变异"):
                    print(f"      父源变异: {result['父源变异']}")
            if result.get("相位证据"):
                print(f"    相位证据句数: {len(result['相位证据'])}")
            if result.get("患者数量"):
                print(f"    患者数量: {result['患者数量']}")

            if result['患者详情']:
                print(f"\n  患者详情:")
                for j, pd in enumerate(result['患者详情'], 1):
                    print(f"    患者{j}: {json.dumps(pd, ensure_ascii=False)}")

            if result.get("变异特征"):
                print(f"\n  变异特征:")
                for k, v in result["变异特征"].items():
                    print(f"    {k}: {v}")

            print(f"\n  功能验证: {result['功能验证']}")
            if result['功能验证详情']:
                for fd in result['功能验证详情']:
                    print(f"    - {fd[:150]}...")

            if result.get("相关句子"):
                print(f"\n  原文相关句:")
                for rs in result["相关句子"][:3]:
                    print(f"    > {rs[:200]}...")

            if result.get("总结段落"):
                print(f"\n  === 标准化总结段落 ===")
                print(f"  {result['总结段落']}")
                print(f"  === 结束 ===")
        else:
            print(f"\n  文献中未提及目标变异 {args.gene} {target_cdna}")
            if result.get("一句话概括"):
                print(f"\n  === 文献简介 ===")
                print(f"  {result['一句话概括']}")
                print(f"  === 结束 ===")

        results.append(result)

        if i < len(unique_pmids) - 1:
            random_sleep()
            request_count += 1
            periodic_long_pause(request_count)

    # 输出JSON
    output_path = os.path.join(os.getcwd(), args.output)
    with open(output_path, "w", encoding="utf-8") as f:
        json.dump(results, f, indent=2, ensure_ascii=False)
    print(f"\n\nJSON output: {output_path}")
    print(f"共处理 {len(results)}/{len(args.pmids)} 篇文献")
    print(f"总请求数: {request_count}")

    # 生成 Excel 汇总表
    excel_dir = args.excel_dir
    safe_gene = re.sub(r'[\\/:*?"<>|]', '_', args.gene) if args.gene else "gene"
    safe_variant = re.sub(r'[\\/:*?"<>|]', '_', target_cdna) if target_cdna else "variant"
    excel_filename = f"{safe_gene}_{safe_variant}_文献汇总.csv"
    excel_path = os.path.join(excel_dir, excel_filename)
    try:
        _generate_excel_csv(results, excel_path)
        print(f"\nExcel 汇总表: {excel_path}")
    except Exception as e:
        print(f"\n[警告] Excel 汇总表生成失败: {e}")

    # 打印最终站点计数
    counts = load_site_counts()
    if counts:
        print(f"\n今日站点访问计数:")
        for site_key, count in counts.items():
            status = "[已达上限]" if count >= DAILY_SITE_LIMIT else "[OK]"
            print(f"   {site_key}: {count}/{DAILY_SITE_LIMIT} {status}")

    if stopped_early:
        print("\n因站点访问限制提前停止")

    # 汇总表格
    print(f"\n{'='*80}")
    print(f"汇总:")
    print(f"{'PMID':<10} {'提及':<6} {'变异类型':<16} {'致病性':<14} {'遗传方式':<14} {'患者数':<6}")
    print("-" * 80)
    for r in results:
        pmid = r.get("PMID", "?")[:10]
        mentioned = "是" if r.get("变异提及") else "否"
        vtype = (r.get("变异类型", "-") or "-")[:16]
        patho = (r.get("致病性", "-") or "-")[:14]
        zyg = (r.get("合子状态", "-") or "-")[:14]
        pcnt = r.get("患者数量", 0) or 0
        print(f"{pmid:<10} {mentioned:<6} {vtype:<16} {patho:<14} {zyg:<14} {pcnt:<6}")


if __name__ == "__main__":
    main()