# PubMed Parser MCP Server

基于 [FastMCP](https://github.com/jlowin/fastmcp) 的 PubMed 文献查询和变异信息提取 MCP 服务器。

## 安装

```bash
cd /Users/roger/Documents/GitHub/mygeno/mcp/pubmed_parser

# 首次安装：创建虚拟环境并安装依赖
uv venv .venv
uv pip install --python .venv/bin/python fastmcp
```

> 注：本项目使用 `uv` 管理虚拟环境，避免 macOS 系统 Python 的 externally-managed 限制。

## 运行

```bash
# 激活虚拟环境
source .venv/bin/activate

# 运行服务器
python server.py
```

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
| `search_pubmed_by_pmid` | 根据 PMID 获取 PubMed 摘要信息 |
| `get_fulltext_by_pmid` | 获取文章全文（优先 PMC 开放获取） |
| `extract_variant_info` | 从文章中提取目标变异信息（致病性、遗传模式、表型等） |
| `search_variant_keywords` | 生成变异搜索关键词变体 |
| `batch_extract_variants` | 批量提取多篇文章的变异信息 |

## 资源

- `pubmed://pmid/{pmid}` — 获取指定 PMID 的文章信息（文本格式）

## 提示 (Prompts)

- `analyze_patient_variants(pmid, gene)` — 生成分析某基因相关患者变异的提示词
