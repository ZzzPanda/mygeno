# FindIt sidecar

Tiny FastAPI service wrapping [metapub.FindIt](https://metapub.readthedocs.io/en/latest/api_findit.html). Resolves a PMID to a publisher PDF URL using metapub's strategies (PMC, Nature, Wiley, Springer, BMC, …). Downloading is the caller's job.

## API

- `GET /healthz` — liveness
- `GET /findit/{pmid}?verify=false&use_nih=true` — returns
  ```json
  {"pmid":"33157158","url":"https://...pdf","reason":null,
   "doi":"10.1038/...","journal":"Nature","source":"nature"}
  ```
  - `url` null when no open PDF (paywall, missing DOI, etc.); `reason` carries metapub's prefix (`PAYWALL:`, `MISSING:`, `TXERROR:`, …).

## Env

| var | default | meaning |
|---|---|---|
| `METAPUB_CACHE_DIR` | `/cache` | metapub SQLite cache (mount a volume to keep it warm) |
| `FINDIT_VERIFY` | `false` | HEAD-check the URL before returning. Slow but accurate. |
| `FINDIT_USE_NIH` | `true` | Prefer NIH/PMC when available |
| `FINDIT_REQUEST_TIMEOUT` | `15` | Per-request timeout passed to metapub |
| `NCBI_API_KEY` | — | (optional) raises NCBI rate limits |
| `LOG_LEVEL` | `INFO` | |

## Run standalone

```bash
docker build -t findit-svc .
docker run --rm -p 8000:8000 -v findit-cache:/cache findit-svc
curl -s http://localhost:8000/findit/33157158 | jq
```

In `pubmed_mcp/docker-compose.yml` it runs as the `findit` service alongside the Rust downloader.
