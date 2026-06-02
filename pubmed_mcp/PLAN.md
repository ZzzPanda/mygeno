# pmid_downloader — Rust MCP Server 实施计划

把 [scihub_downloader.py](../../../../Downloads/skills/scihub-downloader/scripts/scihub_downloader.py)（PMID → DOI → Sci-Hub → 本地 PDF）改造成 MCP 服务，PDF 落 MinIO，元数据落 SQLite，docker-compose 一键起。

---

## 1. 决策一览（已确认）

| 维度 | 选择 | 备注 |
|------|------|------|
| 语言 / 运行时 | Rust + tokio | |
| MCP 框架 | **rmcp**（modelcontextprotocol/rust-sdk 官方）| `rmcp = { version = "*", features = ["server", "transport-streamable-http-server"] }` |
| 传输 | **Streamable HTTP** | 默认监听 `0.0.0.0:8080/mcp`，单端点，适合容器化 |
| 对象存储 | **MinIO**（S3 兼容） | bucket `pmid-pdfs`，key 形如 `PubMed{PMID}.pdf` |
| 元数据库 | **SQLite** | mcp 容器挂卷 `/data/meta.db`，避免再起一个 PG 服务 |
| 工具返回 | **对象存储 URI / 路径** | 形如 `s3://pmid-pdfs/PubMed20981092.pdf` 加元数据；不直接返回 PDF 字节 |
| 部署 | docker-compose | mcp + minio (+ minio-init bucket 引导) |
| 网络 | 直连 Sci-Hub；保留 `HTTP_PROXY` / `HTTPS_PROXY` 环境变量入口 |

**已确认**：

1. 工具集合：`download_pmid` / `download_pmids` / `get_pdf_uri` / `list_downloaded`（`presign_url` 留作后续）。
2. compose 带 `minio-init` 服务做 bucket 自举。

---

## 2. 仓库结构

放在 [pubmed_mcp/](.) 下：

```
pubmed_mcp/
├── Cargo.toml
├── Cargo.lock
├── Dockerfile
├── docker-compose.yml
├── .dockerignore
├── .env.example
├── README.md
├── migrations/
│   └── 0001_init.sql                 # SQLite schema
└── src/
    ├── main.rs                       # 启动 streamable-http server
    ├── config.rs                     # 环境变量 → Config
    ├── http_client.rs                # 共享 reqwest::Client (UA / 代理 / cookie jar)
    ├── pubmed.rs                     # PMID → DOI (E-utilities)
    ├── scihub.rs                     # DOI → PDF bytes（8 级解析 + 完整性校验）
    ├── storage.rs                    # MinIO/S3 put/head/presign
    ├── db.rs                         # SQLite repository (sqlx)
    ├── pipeline.rs                   # 编排：cache hit? → download → store → record
    └── server.rs                     # rmcp ServerHandler + 工具实现
```

---

## 3. 核心数据流

```
download_pmid(pmid)
  ├─► db.lookup(pmid)                       命中 → 直接返回 cached URI
  ├─► pubmed.fetch_doi(pmid)                E-utilities esummary
  ├─► scihub.fetch_pdf(doi)                 镜像轮询 + 8 级 HTML 解析
  ├─► validate_pdf(bytes)                   头 %PDF- + 尾 %%EOF
  ├─► storage.put(key, bytes)               PUT s3://pmid-pdfs/PubMed{pmid}.pdf
  └─► db.insert(pmid, doi, key, size, sha256, downloaded_at)
      返回 { pmid, doi, uri, size_bytes, status: "downloaded" | "cached" }
```

镜像轮询、随机延迟（2–5s 阅读 / 0.8–2.5s 点击 / 3–8s 切镜像 / 6–12s 切 PMID）、CAPTCHA 识别、PDF 8 级链接抽取 —— 全部从 Python 版翻译过来，行为对齐。

---

## 4. SQLite Schema（`migrations/0001_init.sql`）

```sql
CREATE TABLE IF NOT EXISTS downloads (
    pmid          TEXT PRIMARY KEY,
    doi           TEXT NOT NULL,
    object_key    TEXT NOT NULL,            -- 'PubMed{pmid}.pdf'
    bucket        TEXT NOT NULL,
    size_bytes    INTEGER NOT NULL,
    sha256        TEXT NOT NULL,
    source_mirror TEXT NOT NULL,            -- 'sci-hub.ru' / 'st' / 'se'
    downloaded_at INTEGER NOT NULL          -- unix epoch seconds
);
CREATE INDEX IF NOT EXISTS idx_downloads_doi ON downloads(doi);
CREATE INDEX IF NOT EXISTS idx_downloads_at  ON downloads(downloaded_at);
```

---

## 5. MCP 工具定义（rmcp `#[tool]`）

| 工具 | 入参 | 返回 |
|------|------|------|
| `download_pmid` | `pmid: String` | `{ pmid, doi, uri, size_bytes, status }` |
| `download_pmids` | `pmids: Vec<String>`, `concurrency?: u8 = 1` | `Vec<上面的对象>`；遵守人类节奏延迟 |
| `get_pdf_uri` | `pmid: String` | `{ pmid, uri }` 或 `{ pmid, status: "missing" }` |
| `list_downloaded` | `limit?: u32 = 50`, `offset?: u32 = 0` | `Vec<{ pmid, doi, uri, downloaded_at }>` |

`uri` 格式：`s3://{bucket}/{key}`（默认 `s3://pmid-pdfs/PubMed{pmid}.pdf`）。
扩展点：`presign_url(pmid, ttl_secs)` 在二阶段加上，第一阶段不做。

---

## 6. 主要 crate 选型

| 用途 | crate |
|------|-------|
| MCP server | `rmcp` (server + transport-streamable-http-server) |
| 异步运行时 | `tokio` (full) |
| HTTP 客户端 | `reqwest` (rustls, gzip, cookies) |
| HTML 解析 | `scraper`（CSS 选择器，覆盖 8 级策略） |
| 正则（storage 路径） | `regex` |
| S3 SDK | `aws-sdk-s3` + `aws-config`（配 MinIO endpoint / path-style） |
| SQLite | `sqlx`（runtime-tokio-rustls + sqlite + macros） |
| 配置 | `figment` 或裸 `std::env` + `serde` |
| 日志 / trace | `tracing` + `tracing-subscriber` |
| 错误 | `thiserror`（lib 模块） + `anyhow`（main） |
| 哈希 | `sha2` |

---

## 7. 配置（环境变量）

`.env.example`：

```
# MCP
MCP_BIND=0.0.0.0:8080
MCP_PATH=/mcp                       # streamable HTTP 端点

# 存储
S3_ENDPOINT=http://minio:9000
S3_REGION=us-east-1
S3_BUCKET=pmid-pdfs
S3_ACCESS_KEY=minio
S3_SECRET_KEY=minio12345
S3_FORCE_PATH_STYLE=true

# 元数据
SQLITE_PATH=/data/meta.db

# 行为
SCIHUB_MIRRORS=sci-hub.ru,sci-hub.st,sci-hub.se
HTTP_TIMEOUT_SECS=30
DOWNLOAD_TIMEOUT_SECS=120
HUMAN_DELAY=true                    # 复刻 Python 版随机延迟，关掉则全速

# 可选代理
HTTP_PROXY=
HTTPS_PROXY=
```

---

## 8. Dockerfile（多阶段，distroless）

```dockerfile
# --- builder ---
FROM rust:1-bookworm AS builder
WORKDIR /src
COPY Cargo.toml Cargo.lock ./
RUN mkdir src && echo 'fn main(){}' > src/main.rs && cargo build --release && rm -rf src
COPY . .
RUN cargo build --release --bin pmid_downloader

# --- runtime ---
FROM gcr.io/distroless/cc-debian12
COPY --from=builder /src/target/release/pmid_downloader /usr/local/bin/pmid_downloader
COPY --from=builder /src/migrations /migrations
ENV SQLITE_PATH=/data/meta.db
EXPOSE 8080
ENTRYPOINT ["/usr/local/bin/pmid_downloader"]
```

> 如果 distroless 跑 sqlx 迁移有麻烦，回退到 `debian:12-slim` + ca-certificates + 非 root user。

---

## 9. docker-compose.yml

```yaml
services:
  minio:
    image: minio/minio:latest
    command: server /data --console-address ":9001"
    environment:
      MINIO_ROOT_USER: ${S3_ACCESS_KEY:-minio}
      MINIO_ROOT_PASSWORD: ${S3_SECRET_KEY:-minio12345}
    volumes:
      - minio-data:/data
    ports:
      - "9000:9000"      # S3 API
      - "9001:9001"      # web console
    healthcheck:
      test: ["CMD", "curl", "-f", "http://localhost:9000/minio/health/live"]
      interval: 5s
      retries: 20

  minio-init:
    image: minio/mc:latest
    depends_on:
      minio:
        condition: service_healthy
    entrypoint: >
      /bin/sh -c "
      mc alias set local http://minio:9000 ${S3_ACCESS_KEY:-minio} ${S3_SECRET_KEY:-minio12345} &&
      mc mb -p local/${S3_BUCKET:-pmid-pdfs} || true &&
      echo bucket ready
      "

  pmid_downloader:
    build: .
    depends_on:
      minio-init:
        condition: service_completed_successfully
    environment:
      MCP_BIND: 0.0.0.0:8080
      S3_ENDPOINT: http://minio:9000
      S3_BUCKET: ${S3_BUCKET:-pmid-pdfs}
      S3_ACCESS_KEY: ${S3_ACCESS_KEY:-minio}
      S3_SECRET_KEY: ${S3_SECRET_KEY:-minio12345}
      S3_FORCE_PATH_STYLE: "true"
      SQLITE_PATH: /data/meta.db
      HTTP_PROXY: ${HTTP_PROXY:-}
      HTTPS_PROXY: ${HTTPS_PROXY:-}
    volumes:
      - mcp-data:/data
    ports:
      - "8080:8080"

volumes:
  minio-data:
  mcp-data:
```

Claude Code / Desktop 客户端配置：

```json
{
  "mcpServers": {
    "pmid_downloader": {
      "type": "http",
      "url": "http://localhost:8080/mcp"
    }
  }
}
```

---

## 10. 实施阶段（可分别 PR）

1. **骨架** — Cargo 初始化、Dockerfile、compose、health endpoint、空 MCP server 启动。
2. **管道核心** — `pubmed.rs` + `scihub.rs`（含 8 级解析和延迟）+ 行为对齐 Python 版的单元测试。
3. **持久化** — `storage.rs`（MinIO put/head）+ `db.rs`（sqlx 迁移 + repo）+ `pipeline.rs`。
4. **MCP 暴露** — `server.rs` 注册四个工具，端到端联调（mc 看 bucket、`sqlite3 meta.db` 看记录）。
5. **打磨** — 错误码标准化、tracing 日志、README、`.env.example`、CI（fmt / clippy / test / docker build）。

---

## 11. 风险点 & 应对

| 风险 | 应对 |
|------|------|
| Sci-Hub 镜像被墙或返 CAPTCHA | 保留 `HTTP_PROXY` 入口；CAPTCHA 检测沿用 Python 版（`recaptcha` / `cf-captcha` / `cf-wrapper`） |
| `aws-sdk-s3` 体积大、编译慢 | 二选一：(a) 接受首次构建慢，多阶段 + cargo-chef 缓存依赖；(b) 换 `rust-s3` 更轻 —— 若构建时间不可接受再切 |
| sqlx 编译期校验 vs 容器构建 | 用 `query!` 时需要 `DATABASE_URL` 或离线数据；首版可用 `query` (运行时校验) 避开离线模式 |
| PDF 命名碰撞 / 大小写差异 | 一律用 `PubMed{pmid}.pdf`，`pmid` 已 `.is_digit()` 校验 |
| MinIO 密钥默认值入仓 | `.env.example` 给样例，真值从 `.env`（gitignore）注入 |

---

## 12. Definition of Done

- [ ] `docker compose up --build` 能起来，3 个 service 全部 healthy。
- [ ] Claude Desktop 接入 `http://localhost:8080/mcp`，能列出 4 个工具。
- [ ] `download_pmid` 跑通一个公开样例（如 `20981092`），MinIO 看到对象、SQLite 看到记录、第二次调用走缓存。
- [ ] `download_pmids` 多个 PMID，行为节奏符合 Python 版（随机延迟、镜像轮询）。
- [ ] README 含本地启动 / 客户端配置 / 故障排查（CAPTCHA、proxy、bucket 不存在）。
