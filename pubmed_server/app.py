"""FindIt sidecar — wraps metapub.FindIt behind a small HTTP API.

Resolves a PMID to a publisher PDF URL via metapub's "dance" strategies
(open-access journals, PMC, Nature, Wiley, etc.). Returns the URL plus
provenance; the caller is responsible for downloading the PDF.
"""

from __future__ import annotations

import logging
import os
from typing import Optional

from fastapi import FastAPI, HTTPException, Query
from pydantic import BaseModel

from metapub import FindIt
from metapub.exceptions import AccessDenied, NoPDFLink

logging.basicConfig(level=os.getenv("LOG_LEVEL", "INFO"))
logger = logging.getLogger("findit-svc")

CACHE_DIR = os.getenv("METAPUB_CACHE_DIR", "/cache")
DEFAULT_VERIFY = os.getenv("FINDIT_VERIFY", "false").lower() in ("true", "1", "yes")
DEFAULT_USE_NIH = os.getenv("FINDIT_USE_NIH", "true").lower() in ("true", "1", "yes")
REQUEST_TIMEOUT = int(os.getenv("FINDIT_REQUEST_TIMEOUT", "15"))

os.makedirs(CACHE_DIR, exist_ok=True)

app = FastAPI(title="FindIt sidecar", version="0.1.0")


class FindItResult(BaseModel):
    pmid: str
    url: Optional[str] = None
    reason: Optional[str] = None
    doi: Optional[str] = None
    journal: Optional[str] = None
    source: Optional[str] = None


def _classify(url: str) -> str:
    u = url.lower()
    for needle, label in [
        ("ncbi.nlm.nih.gov/pmc", "pmc"),
        ("europepmc.org", "epmc"),
        ("nature.com", "nature"),
        ("sciencedirect", "sciencedirect"),
        ("wiley.com", "wiley"),
        ("springer", "springer"),
        ("biomedcentral", "bmc"),
        ("plos.org", "plos"),
        ("bmj.com", "bmj"),
        ("oup.com", "oup"),
        ("sagepub", "sage"),
        ("jamanetwork", "jama"),
        ("cell.com", "cell"),
        ("nih.gov", "nih"),
    ]:
        if needle in u:
            return label
    return "other"


@app.get("/healthz")
def healthz() -> dict:
    return {"ok": True}


@app.get("/findit/{pmid}", response_model=FindItResult)
def findit(
    pmid: str,
    verify: Optional[bool] = Query(default=None),
    use_nih: Optional[bool] = Query(default=None),
) -> FindItResult:
    if not pmid.isdigit():
        raise HTTPException(status_code=400, detail="pmid must be all digits")

    v = DEFAULT_VERIFY if verify is None else verify
    n = DEFAULT_USE_NIH if use_nih is None else use_nih

    try:
        src = FindIt(
            pmid,
            verify=v,
            use_nih=n,
            cachedir=CACHE_DIR,
            request_timeout=REQUEST_TIMEOUT,
        )
    except (NoPDFLink, AccessDenied) as e:
        return FindItResult(pmid=pmid, url=None, reason=str(e))
    except Exception as e:
        logger.exception("findit failed for pmid=%s", pmid)
        raise HTTPException(status_code=502, detail=f"findit error: {e}") from e

    pma = src.pma
    doi = getattr(pma, "doi", None) if pma is not None else None
    journal = getattr(pma, "journal", None) if pma is not None else None

    return FindItResult(
        pmid=pmid,
        url=src.url,
        reason=src.reason,
        doi=doi,
        journal=journal,
        source=_classify(src.url) if src.url else None,
    )
