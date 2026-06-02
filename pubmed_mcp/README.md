# pmid_downloader

Rust MCP server that downloads PubMed full-text PDFs via Sci-Hub. PDFs land in MinIO (S3-compatible), metadata in SQLite. Tools return `s3://` URIs.

Pipeline: `PMID → DOI (PubMed E-utils) → Sci-Hub mirror rotation → 8-strategy HTML parse → %PDF- / %%EOF validation → S3 PUT → SQLite insert`.

## Quick start

```bash
cp .env.example .env
docker compose up --build
```

Three services come up: `minio` (S3 API on `:9000`, console on `:9001`), `minio-init` (one-shot bucket creator), and `pmid_downloader` (MCP on `:8080/mcp`).

Add to your MCP client config:

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

## Tools

| Tool | Args | Returns |
|------|------|---------|
| `download_pmid` | `pmid: string` | `{pmid, doi, uri, size_bytes, status}` — `cached` on repeat |
| `download_pmids` | `pmids: string[]` | `[{...}]` — sequential with human-paced delays |
| `get_pdf_uri` | `pmid: string` | `{pmid, uri, status}` — lookup only, no fetch |
| `list_downloaded` | `limit?, offset?` | `[{pmid, doi, uri, size_bytes, downloaded_at}]` |

`uri` format: `s3://pmid-pdfs/PubMed{pmid}.pdf`.

## Configuration

All via env vars; see [`.env.example`](.env.example). Highlights:

- `SCIHUB_MIRRORS` — comma-separated, tried in order
- `HUMAN_DELAY` — set `false` to disable random pacing (faster, riskier)
- `HTTP_PROXY` / `HTTPS_PROXY` — optional egress proxy (some regions can't reach Sci-Hub directly)

## Local dev

```bash
cargo run
```

Requires a reachable MinIO. Easiest: `docker compose up minio minio-init -d`, then run the binary on the host:

```bash
S3_ENDPOINT=http://localhost:9000 SQLITE_PATH=./meta.db cargo run
```

## Troubleshooting

- **`captcha/blocked` in logs** → mirror returned a CAPTCHA page; loop tries the next. Configure a proxy if all mirrors fail.
- **`no_doi` status** → PubMed has no DOI for that PMID; nothing to do.
- **`not_found`** → DOI exists but Sci-Hub doesn't have the article (often very new papers).
- **Bucket missing** → the binary auto-creates it on startup; `minio-init` is a belt-and-suspenders pre-seed.
