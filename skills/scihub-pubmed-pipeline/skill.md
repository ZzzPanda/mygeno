---
name: scihub-pubmed-pipeline
description: 全流程管道：Sci-Hub 下载 PDF → 离线 PDF 变异搜索 → 在线 PubMed API 提取。输入 PMID + 基因 + 变异，一键完成从下载到解读的全流程。输出 JSON + Excel 汇总表至 D:\claude_code\project1\文献提取结果。
---

# Sci-Hub → PubMed 全流程管道

一键完成三步骤：**Sci-Hub 下载 PDF** → **离线 PDF 变异搜索** → **在线 PubMed API 提取**。

## 工作流程

| 步骤 | 操作 | 脚本 | 说明 |
|------|------|------|------|
| 1 | Sci-Hub 下载 | [scihub_downloader.py](../scihub-downloader/scripts/scihub_downloader.py) | PMID → DOI → Sci-Hub → PDF |
| 2 | 离线 PDF 搜索 | [pdf_variant_search.py](../pubmed-extractor/scripts/pdf_variant_search.py) | 本地 PDF 全文变异搜索 + 表型/致病性/相位提取 |
| 3 | 在线 API 提取 | [pubmed_extractor.py](../pubmed-extractor/scripts/pubmed_extractor.py) | NCBI + Europe PMC + PMC 全文本 API 提取 |

## 使用方法

### 环境准备（首次使用）

```bash
pip install requests beautifulsoup4 pdfplumber
```

### 一行命令执行

```bash
# 完整流程（下载 + PDF 搜索 + 在线提取）
python .claude/skills/scihub-pubmed-pipeline/scripts/scihub_pubmed_pipeline.py \
    --pmids 20981092 --gene ABCA4 --variant "c.763C>T (p.Arg255Cys)"

# 多 PMID + 转录本
python .claude/skills/scihub-pubmed-pipeline/scripts/scihub_pubmed_pipeline.py \
    --pmids 20981092 33090715 33301772 \
    --gene PKD1 \
    --variant "c.1522T>C (p.Cys508Arg)" \
    --transcript NM_001009944.3

# 跳过下载（PDF 已存在，仅做搜索 + 在线提取）
python .claude/skills/scihub-pubmed-pipeline/scripts/scihub_pubmed_pipeline.py \
    --pmids 20981092 --gene ABCA4 --variant "c.763C>T (p.Arg255Cys)" \
    --skip-download

# 仅下载 + PDF 搜索（跳过在线 API，节省时间）
python .claude/skills/scihub-pubmed-pipeline/scripts/scihub_pubmed_pipeline.py \
    --pmids 20981092 --gene ABCA4 --variant "c.763C>T (p.Arg255Cys)" \
    --skip-online
```

### 参数说明

| 参数 | 必需 | 说明 |
|------|------|------|
| `--pmids` | 是 | PubMed ID 列表，空格分隔 |
| `--gene` | 是 | 目标基因名称 |
| `--variant` | 是 | 目标变异，格式如 `"c.763C>T (p.Arg255Cys)"` |
| `--transcript` | 否 | 转录本 ID（推荐提供，可识别版本差异） |
| `--pdf-dir` | 否 | PDF 下载/搜索目录，默认 `D:\claude_code\project1\sci` |
| `--excel-dir` | 否 | Excel 汇总表输出目录，默认 `D:\claude_code\project1\文献提取结果` |
| `--skip-download` | 否 | 跳过 Sci-Hub 下载步骤（PDF 已存在时使用） |
| `--skip-online` | 否 | 跳过在线 PubMed API 提取步骤（仅需 PDF 搜索时使用） |

## 输出文件

### 步骤 1 — PDF 下载
- 路径: `D:\claude_code\project1\sci\PubMed{PMID}.pdf`
- 已存在的文件自动跳过

### 步骤 2 — 离线 PDF 搜索结果
- JSON: `{excel-dir}\{基因}_pdf_search_results.json`
- Excel: `{excel-dir}\{基因}_{变异}_PDF文献汇总.csv`

| 列名 | 说明 |
|------|------|
| PMID | PubMed ID |
| 标题 | PDF 前 300 字符提取的标题 |
| 是否提及此位点 | 是 / 否 |
| 患者数 | 估算的携带目标变异患者数 |
| 致病性 | 致病性评级 |
| 关联合子状态 | 纯合/杂合/复合杂合 |
| 临床表型 | 提取的表型关键词 |
| 匹配关键词 | 实际匹配到的变异关键词 |
| 相关句子片段 | PDF 中包含变异提及的句子 |

### 步骤 3 — 在线 API 提取结果
- JSON: `{excel-dir}\{基因}_online_results.json`
- Excel: `{excel-dir}\{基因}_{变异}_文献汇总.csv`

详细字段参见 [pubmed-extractor SKILL.md](../pubmed-extractor/SKILL.md#输出文件)。

## 各步骤独立性

每个步骤可独立跳过，适应不同场景：

| 场景 | 命令 |
|------|------|
| 全新文献，需下载 + 全面分析 | 完整流程（无 skip 参数） |
| PDF 已下载过，仅需分析 | `--skip-download` |
| 仅需 PDF 内容搜索，不需在线 API | `--skip-online` |
| 仅需在线 API 提取 | 直接使用 [pubmed-extractor](../pubmed-extractor/) |

## 限制

- Sci-Hub 未收录的文章无法下载（特别是较新的文章）
- 需要网络能访问 Sci-Hub 镜像站和 PubMed API
- PDF 搜索依赖 pdfplumber 提取文本质量（扫描版 PDF 效果差）
- 在线 API 对非开放获取文献仅能获取摘要
- 各步骤速率限制独立运作，完整流程耗时较长

## 依赖脚本

| 脚本 | 路径 |
|------|------|
| 管道主脚本 | `.claude/skills/scihub-pubmed-pipeline/scripts/scihub_pubmed_pipeline.py` |
| Sci-Hub 下载器 | `.claude/skills/scihub-downloader/scripts/scihub_downloader.py` |
| PDF 变异搜索 | `.claude/skills/pubmed-extractor/scripts/pdf_variant_search.py` |
| PubMed 在线提取 | `.claude/skills/pubmed-extractor/scripts/pubmed_extractor.py` |