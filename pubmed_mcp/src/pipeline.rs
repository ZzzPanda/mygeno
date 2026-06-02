use std::{
    sync::Arc,
    time::{SystemTime, UNIX_EPOCH},
};

use rand::Rng;
use reqwest::Client;
use sha2::{Digest, Sha256};
use tracing::{info, warn};

use crate::{
    config::Config,
    db::{Db, DownloadRow},
    pubmed, scihub,
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

        let doi = match pubmed::fetch_doi(&self.http, &pmid).await {
            Ok(Some(d)) => d,
            Ok(None) => {
                return DownloadOutcome {
                    pmid,
                    doi: None,
                    uri: None,
                    size_bytes: None,
                    status: "no_doi".into(),
                    error: Some("PubMed has no DOI for this PMID".into()),
                }
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
        };

        let hit = match scihub::fetch_pdf(&self.http, &self.cfg, &doi).await {
            Ok(h) => h,
            Err(e) => {
                return DownloadOutcome {
                    pmid,
                    doi: Some(doi),
                    uri: None,
                    size_bytes: None,
                    status: "not_found".into(),
                    error: Some(e.to_string()),
                };
            }
        };

        if let Err(msg) = scihub::validate_pdf(&hit.bytes) {
            return DownloadOutcome {
                pmid,
                doi: Some(doi),
                uri: None,
                size_bytes: None,
                status: "invalid_pdf".into(),
                error: Some(msg.into()),
            };
        }

        let key = format!("PubMed{pmid}.pdf");
        if let Err(e) = self.storage.put(&key, hit.bytes.clone()).await {
            return DownloadOutcome {
                pmid,
                doi: Some(doi),
                uri: None,
                size_bytes: None,
                status: "error".into(),
                error: Some(format!("S3 put: {e}")),
            };
        }

        let mut hasher = Sha256::new();
        hasher.update(&hit.bytes);
        let sha256 = hex::encode(hasher.finalize());

        let row = DownloadRow {
            pmid: pmid.clone(),
            doi: doi.clone(),
            object_key: key.clone(),
            bucket: self.storage.bucket().to_string(),
            size_bytes: hit.bytes.len() as i64,
            sha256,
            source_mirror: hit.mirror,
            downloaded_at: SystemTime::now()
                .duration_since(UNIX_EPOCH)
                .map(|d| d.as_secs() as i64)
                .unwrap_or(0),
        };
        if let Err(e) = self.db.insert(&row).await {
            warn!(error = %e, "db insert failed (object already in S3)");
        }

        info!(%pmid, %doi, size = row.size_bytes, "downloaded");
        DownloadOutcome {
            pmid,
            doi: Some(doi),
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
