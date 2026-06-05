# PubMed Parser MCP Server

基于 [FastMCP](https://github.com/jlowin/fastmcp) 的 PDF 文献解析和变异信息提取 MCP 服务器。

## 安装

```bash
cd /Users/roger/Documents/GitHub/mygeno/mcp/pubmed_parser

# 创建虚拟环境并安装依赖
uv venv .venv
uv pip install fastmcp pymupdf pytest
```

## 运行

### 测试模式（验证安装）
```bash
.venv/bin/python server.py --test
```

### 启动 MCP 服务器
MCP 服务器通过 stdio 通信，需要通过 Claude Code MCP 配置连接：

```bash
.venv/bin/python server.py
```

直接运行会显示 JSON 解析错误，这是正常现象（没有 MCP 客户端连接）。

## Claude Code 集成

在 Claude Code 的 MCP 配置中添加：

```json
{
  "mcpServers": {
    "pubmed-parser": {
      "command": "/Users/roger/Documents/GitHub/mygeno/mcp/pubmed_parser/.venv/bin/python",
      "args": ["/Users/roger/Documents/GitHub/mygeno/mcp/pubmed_parser/server.py"]
    }
  }
}
```

## 提供的工具

| 工具 | 说明 |
|------|------|
| `extract_variant_info` | 从 PDF 文件中提取目标变异信息（致病性、遗传模式、表型等） |
| `parse_pdf` | 解析 PDF 文件，提取文本和表格 |
| `analyze_variant` | 从文本/摘要中提取变异信息（无需 PDF） |
| `search_variant_keywords` | 生成变异搜索关键词的所有变体 |

## 使用示例

```python
# 从 PDF 提取变异信息
extract_variant_info(
    pdf_path="/path/to/PMID:33374015.pdf",
    cdna="c.2279C>T",
    protein="p.Thr760Met",
    gene="CFTR",
    transcript="NM_000492.4"
)

# 解析 PDF 获取文本和表格
parse_pdf("/path/to/PMID:33374015.pdf")

# 直接从文本分析变异
analyze_variant(
    text="The patient carried the c.1166G>A variant...",
    cdna="c.1166G>A",
    protein="p.R389H",
    gene="GENE"
)
```

## 资源

- `pdf://parse/{pdf_path}` — 解析指定路径的 PDF 并返回文本内容

## 提示 (Prompts)

- `analyze_patient_variants(pdf_path, gene)` — 生成分析某基因相关患者变异的提示词

## 运行测试

```bash
# pytest（所有测试）
.venv/bin/python -m pytest tests/test_pdf_parser.py -v

# 单个 PDF 测试（调试用）
.venv/bin/python -m pytest tests/test_pdf_parser.py::TestPDFParser::test_single_pdf -v

# CLI 模式
.venv/bin/python tests/test_pdf_parser.py --data-dir data/pdf
```
