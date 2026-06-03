---
name: scihub-downloader
description: 根据 PMID 从 Sci-Hub 下载 PDF 全文。流程：PMID → DOI (PubMed API) → Sci-Hub → PDF。依次尝试 sci-hub.ru → sci-hub.st → sci-hub.se 三个镜像站，内置随机延迟模拟人类操作节奏。保存至 D:\claude_code\project1\sci，命名为 PubMed{PMID}.pdf。
---

# Sci-Hub PDF 下载器

根据 PubMed ID (PMID) 通过 Sci-Hub 自动下载 PDF 全文。

## 工作流程

| 步骤 | 操作 | 说明 |
|------|------|------|
| 1 | PMID → DOI | 通过 PubMed E-utilities API 获取 DOI |
| 2 | Sci-Hub 访问 | `http://{mirror}/{DOI}` 访问文章页 |
| 3 | PDF 链接提取 | 8 级策略识别嵌入的 PDF 链接 |
| 4 | PDF 下载 | 带 Referer + 随机延迟下载 |
| 5 | 完整性验证 | 检查 `%PDF-` 文件头和 `%%EOF` 尾部 |

## 使用方法

### 环境准备（首次使用）

```bash
pip install requests beautifulsoup4
```

### 下载

```bash
# 单个 PMID
python .claude/skills/scihub-downloader/scripts/scihub_downloader.py 20981092

# 多个 PMID
python .claude/skills/scihub-downloader/scripts/scihub_downloader.py 20981092 33090715 33301772

# 逗号分隔
python .claude/skills/scihub-downloader/scripts/scihub_downloader.py --pmids 20981092,33090715

# 自定义输出目录
python .claude/skills/scihub-downloader/scripts/scihub_downloader.py 20981092 --output D:\my_pdfs
```

## 访问策略（模拟人类操作）

所有延迟均为随机小数，禁止固定整数：

| 场景 | 延迟范围 | 说明 |
|------|---------|------|
| 打开页面 → 操作 | 2.0 ~ 5.0s | 模拟阅读/反应时间 |
| 看到按钮 → 点击 | 0.8 ~ 2.5s | 模拟点击反应 |
| 连续请求间隔 | 3.0 ~ 8.0s | 模拟操作间隔 |
| 一轮结束 → 下一轮 | 6.0 ~ 12.0s | 模拟休息时间 |

启动时自动访问 Sci-Hub 首页获取会话 cookie（预热），后续请求复用同一会话。

## 镜像站优先级

| 优先级 | 镜像 | 说明 |
|--------|------|------|
| 1 | sci-hub.ru | 主要镜像，HTTP/HTTPS 双协议 |
| 2 | sci-hub.st | 备用镜像 |
| 3 | sci-hub.se | 备用镜像 |

## PDF 链接识别策略（8 级回退）

| 优先级 | 策略 | 说明 |
|--------|------|------|
| 1 | `<meta name="citation_pdf_url">` | 文献标准元数据 |
| 2 | `<object type="application/pdf">` | PDF 嵌入对象 |
| 3 | `<div class="download">` 内链接 | Sci-Hub 下载按钮 |
| 4 | `<iframe src="...">` | 内嵌框架 |
| 5 | `<embed src="...">` | 嵌入式内容 |
| 6 | `<button onclick="...">` | 按钮点击跳转 |
| 7 | `/storage/.../*.pdf` 正则 | Sci-Hub 存储路径 |
| 8 | 任意 `.pdf` 链接 | 通用回退 |

## 输出

- 命名规则: `PubMed{PMID}.pdf`
- 输出目录: `D:\claude_code\project1\sci`
- 已存在的文件自动跳过

## 限制

- Sci-Hub 未收录的文章无法下载（特别是较新的文章）
- 需要网络能访问 Sci-Hub 镜像站
- 频繁请求可能触发限流（脚本已内置随机延迟缓解）