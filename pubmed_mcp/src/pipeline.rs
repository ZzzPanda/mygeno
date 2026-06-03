use std::{
    sync::Arc,
    time::{SystemTime, UNIX_EPOCH},
};

use rand::Rng;
use reqwest::Client;
use sha2::{Digest, Sha256};
use tracing::{debug, info, warn};

use crate::{
    config::Config,
    db::{Db, DownloadRow},
    findit, pmc, pubmed, scihub,
    storage::Storage,
};

#[derive(Debug, Clone, serde::Serialize, schemars::JsonSchema)]
pub struct DownloadOutcome {
    pub pmid: String,
    pub doi: Option<String>,
    pub uri: Option<String>,
    pub size_bytes: Option<i64>,
    pub status: String,
    pub error: Option<String>,
}

#[derive(Clone)]
pub struct Pipeline {
    pub cfg: Config,
    pub http: Arc<Client>,
    pub db: Db,
    pub storage: Storage,
}

impl Pipeline {
    pub async fn download_one(&self, pmid: &str) -> DownloadOutcome {
        let pmid = pmid.trim().to_string();
        if pmid.is_empty() || !pmid.chars().all(|c| c.is_ascii_digit()) {
            return DownloadOutcome {
                pmid,
                doi: None,
                uri: None,
                size_bytes: None,
                status: "invalid".into(),
                error: Some("PMID must be all digits".into()),
            };
        }

        if let Ok(Some(row)) = self.db.lookup(&pmid).await {
            return DownloadOutcome {
                pmid: row.pmid,
                doi: Some(row.doi),
                uri: Some(self.storage.uri_for(&row.object_key)),
                size_bytes: Some(row.size_bytes),
                status: "cached".into(),
                error: None,
            };
        }

        // Try FindIt first (open-access publishers + PMC). Captures DOI as a side effect.
        let mut doi: Option<String> = None;
        let mut findit_err: Option<String> = None;
        let findit_hit = if self.cfg.findit_url.is_some() {
            match findit::fetch_pdf(&self.http, &self.cfg, &pmid).await {
                Ok(hit) => {
                    info!(%pmid, source = %hit.source, "got PDF from FindIt");
                    if let Some(d) = &hit.doi {
                        doi = Some(d.clone());
                    }
                    Some((hit.bytes, hit.source))
                }
                Err(e) => {
                    debug!(%pmid, error = %e, "FindIt unavailable, falling back");
                    findit_err = Some(e.to_string());
                    None
                }
            }
        } else {
            None
        };

        let (bytes, source) = if let Some(hit) = findit_hit {
            hit
        } else {
            match pmc::fetch_pdf(&self.http, &self.cfg, &pmid).await {
                Ok(hit) => {
                    info!(%pmid, source = %hit.source, "got PDF from Europe PMC");
                    (hit.bytes, hit.source)
                }
                Err(e) => {
                    debug!(%pmid, error = %e, "PMC unavailable, falling back to Sci-Hub");
                    if doi.is_none() {
                        match pubmed::fetch_doi(&self.http, &pmid).await {
                            Ok(Some(d)) => doi = Some(d),
                            Ok(None) => {
                                let mut err = String::from("no DOI; ");
                                if let Some(fe) = &findit_err {
                                    err.push_str(&format!("findit: {fe}; "));
                                }
                                err.push_str(&format!("pmc: {e}"));
                                return DownloadOutcome {
                                    pmid,
                                    doi: None,
                                    uri: None,
                                    size_bytes: None,
                                    status: "no_doi".into(),
                                    error: Some(err),
                                };
                            }
                            Err(e) => {
                                warn!(error = %e, "pubmed query failed");
                                return DownloadOutcome {
                                    pmid,
                                    doi: None,
                                    uri: None,
                                    size_bytes: None,
                                    status: "error".into(),
                                    error: Some(format!("PubMed: {e}")),
                                };
                            }
                        }
                    }
                    let doi_val = doi.clone().unwrap();
                    match scihub::fetch_pdf(&self.http, &self.cfg, &doi_val).await {
                        Ok(hit) => (hit.bytes, format!("scihub:{}", hit.mirror)),
                        Err(sh_err) => {
                            let mut err = String::new();
                            if let Some(fe) = &findit_err {
                                err.push_str(&format!("findit: {fe}; "));
                            }
                            err.push_str(&format!("pmc: {e}; scihub: {sh_err}"));
                            return DownloadOutcome {
                                pmid,
                                doi: Some(doi_val),
                                uri: None,
                                size_bytes: None,
                                status: "not_found".into(),
                                error: Some(err),
                            };
                        }
                    }
                }
            }
        };

        // Backfill DOI if Sci-Hub path wasn't taken and FindIt didn't supply one.
        if doi.is_none() {
            if let Ok(Some(d)) = pubmed::fetch_doi(&self.http, &pmid).await {
                doi = Some(d);
            }
        }

        if let Err(msg) = scihub::validate_pdf(&bytes) {
            return DownloadOutcome {
                pmid,
                doi,
                uri: None,
                size_bytes: None,
                status: "invalid_pdf".into(),
                error: Some(msg.into()),
            };
        }

        let key = format!("PubMed{pmid}.pdf");
        if let Err(e) = self.storage.put(&key, bytes.clone()).await {
            return DownloadOutcome {
                pmid,
                doi,
                uri: None,
                size_bytes: None,
                status: "error".into(),
                error: Some(format!("S3 put: {e}")),
            };
        }

        let mut hasher = Sha256::new();
        hasher.update(&bytes);
        let sha256 = hex::encode(hasher.finalize());

        let row = DownloadRow {
            pmid: pmid.clone(),
            doi: doi.clone().unwrap_or_default(),
            object_key: key.clone(),
            bucket: self.storage.bucket().to_string(),
            size_bytes: bytes.len() as i64,
            sha256,
            source_mirror: source,
            downloaded_at: SystemTime::now()
                .duration_since(UNIX_EPOCH)
                .map(|d| d.as_secs() as i64)
                .unwrap_or(0),
        };
        if let Err(e) = self.db.insert(&row).await {
            warn!(error = %e, "db insert failed (object already in S3)");
        }

        info!(%pmid, doi = ?doi, size = row.size_bytes, source = %row.source_mirror, "downloaded");
        DownloadOutcome {
            pmid,
            doi,
            uri: Some(self.storage.uri_for(&key)),
            size_bytes: Some(row.size_bytes),
            status: "downloaded".into(),
            error: None,
        }
    }

    pub async fn download_many(&self, pmids: &[String]) -> Vec<DownloadOutcome> {
        let mut out = Vec::with_capacity(pmids.len());
        for (i, pmid) in pmids.iter().enumerate() {
            if i > 0 && self.cfg.human_delay {
                let secs = rand::thread_rng().gen_range(6.0..12.0);
                tokio::time::sleep(std::time::Duration::from_secs_f64(secs)).await;
            }
            out.push(self.download_one(pmid).await);
        }
        out
    }
}
