#!/usr/bin/env python3
"""
Sci-Hub -> PubMed 全流程管道：下载 PDF -> 离线 PDF 搜索 -> 在线 API 提取。
三步串联，一键完成。

用法:
  python scihub_pubmed_pipeline.py \
      --pmids 20981092 33090715 \
      --gene ABCA4 \
      --variant "c.763C>T (p.Arg255Cys)" \
      [--transcript NM_000350.3] \
      [--pdf-dir D:/claude_code/project1/sci] \
      [--excel-dir D:/claude_code/project1/文献提取结果] \
      [--skip-download] \
      [--skip-online]
"""

import argparse
import os
import subprocess
import sys

# Windows UTF-8
if sys.platform == 'win32':
    sys.stdout.reconfigure(encoding='utf-8')
    sys.stderr.reconfigure(encoding='utf-8')

# 项目根目录和脚本路径
PROJECT_ROOT = r"D:\claude_code\project1"
SKILLS_DIR = os.path.join(PROJECT_ROOT, ".claude", "skills")

SCIHUB_SCRIPT = os.path.join(SKILLS_DIR, "scihub-downloader", "scripts", "scihub_downloader.py")
PDF_SEARCH_SCRIPT = os.path.join(SKILLS_DIR, "pubmed-extractor", "scripts", "pdf_variant_search.py")
PUBMED_EXTRACTOR_SCRIPT = os.path.join(SKILLS_DIR, "pubmed-extractor", "scripts", "pubmed_extractor.py")

DEFAULT_PDF_DIR = os.path.join(PROJECT_ROOT, "sci")
DEFAULT_EXCEL_DIR = os.path.join(PROJECT_ROOT, "文献提取结果")


def run_step(step_name, cmd):
    """运行一个子进程步骤，失败时返回 False。"""
    print(f"\n{'='*60}")
    print(f">>> 步骤: {step_name}")
    print(f">>> 命令: {' '.join(cmd)}")
    print(f"{'='*60}\n")
    result = subprocess.run(cmd, cwd=PROJECT_ROOT)
    if result.returncode != 0:
        print(f"\n!!! 步骤失败 (exit code {result.returncode}): {step_name}")
        return False
    print(f"\n✓ 步骤完成: {step_name}")
    return True


def main():
    parser = argparse.ArgumentParser(
        description="Sci-Hub 下载 → PDF 搜索 → PubMed 在线提取 全流程管道",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
示例:
  # 完整流程（下载 + PDF 搜索 + 在线提取）
  python scihub_pubmed_pipeline.py --pmids 20981092 --gene ABCA4 --variant "c.763C>T (p.Arg255Cys)"

  # 多 PMID
  python scihub_pubmed_pipeline.py --pmids 20981092 33090715 33301772 --gene PKD1 --variant "c.1522T>C (p.Cys508Arg)" --transcript NM_001009944.3

  # 跳过下载（PDF 已存在）
  python scihub_pubmed_pipeline.py --pmids 20981092 --gene ABCA4 --variant "c.763C>T (p.Arg255Cys)" --skip-download

  # 仅下载 + PDF 搜索（跳过在线 API）
  python scihub_pubmed_pipeline.py --pmids 20981092 --gene ABCA4 --variant "c.763C>T (p.Arg255Cys)" --skip-online
        """,
    )
    parser.add_argument("--pmids", nargs="+", required=True, help="PubMed ID 列表，空格分隔")
    parser.add_argument("--gene", required=True, help="目标基因名称")
    parser.add_argument("--variant", required=True, help='目标变异，格式如 "c.763C>T (p.Arg255Cys)"')
    parser.add_argument("--transcript", default="", help="转录本 ID（推荐提供）")
    parser.add_argument("--pdf-dir", default=DEFAULT_PDF_DIR, help=f"PDF 下载/搜索目录 (默认: {DEFAULT_PDF_DIR})")
    parser.add_argument("--excel-dir", default=DEFAULT_EXCEL_DIR, help=f"Excel 输出目录 (默认: {DEFAULT_EXCEL_DIR})")
    parser.add_argument("--skip-download", action="store_true", help="跳过 Sci-Hub 下载步骤")
    parser.add_argument("--skip-online", action="store_true", help="跳过在线 PubMed API 提取步骤")

    args = parser.parse_args()
    pmids = args.pmids

    print(f"Sci-Hub → PubMed 全流程管道")
    print(f"PMID: {', '.join(pmids)}")
    print(f"基因: {args.gene}")
    print(f"变异: {args.variant}")
    if args.transcript:
        print(f"转录本: {args.transcript}")
    print(f"PDF 目录: {args.pdf_dir}")
    print(f"Excel 目录: {args.excel_dir}")

    # ---- 步骤 1: Sci-Hub 下载 ----
    if not args.skip_download:
        download_cmd = [
            sys.executable, SCIHUB_SCRIPT,
            "--pmids", ','.join(pmids),
            "--output", args.pdf_dir,
        ]
        if not run_step("1/3 Sci-Hub PDF 下载", download_cmd):
            print("下载步骤失败，将继续尝试搜索已有 PDF...")
    else:
        print("\n[跳过] 步骤 1/3: Sci-Hub PDF 下载")

    # ---- 步骤 2: 离线 PDF 变异搜索 ----
    pdf_output = os.path.join(args.excel_dir, f"{args.gene}_pdf_search_results.json")
    pdf_search_cmd = [
        sys.executable, PDF_SEARCH_SCRIPT,
        "--pdf-dir", args.pdf_dir,
        "--gene", args.gene,
        "--variant", args.variant,
        "--output", pdf_output,
        "--excel-dir", args.excel_dir,
    ]
    if args.transcript:
        pdf_search_cmd.extend(["--transcript", args.transcript])

    if not run_step("2/3 离线 PDF 变异搜索", pdf_search_cmd):
        print("PDF 搜索步骤失败")

    # ---- 步骤 3: 在线 PubMed API 提取 ----
    if not args.skip_online:
        online_output = os.path.join(args.excel_dir, f"{args.gene}_online_results.json")
        online_cmd = [
            sys.executable, PUBMED_EXTRACTOR_SCRIPT,
            "--pmids"] + pmids + [
            "--gene", args.gene,
            "--variant", args.variant,
            "--output", online_output,
            "--excel-dir", args.excel_dir,
        ]
        if args.transcript:
            online_cmd.extend(["--transcript", args.transcript])

        # v9: 传入PDF搜索结果用于交叉引用（修复摘要级文献遗漏变异检出）
        online_cmd.extend(["--pdf-results", pdf_output])

        if not run_step("3/3 在线 PubMed API 提取", online_cmd):
            print("在线提取步骤失败")
    else:
        print("\n[跳过] 步骤 3/3: 在线 PubMed API 提取")

    # ---- 完成 ----
    print(f"\n{'='*60}")
    print("全流程完成!")
    print(f"  PDF 下载目录: {args.pdf_dir}")
    print(f"  PDF 搜索结果: {pdf_output}")
    print(f"  Excel 汇总表:  {args.excel_dir}")
    if not args.skip_online:
        print(f"  在线提取结果: {online_output}")
    print(f"{'='*60}")


if __name__ == "__main__":
    main()