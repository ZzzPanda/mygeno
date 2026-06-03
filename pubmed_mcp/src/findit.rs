use std::sync::Arc;

use reqwest::Client;
use serde::Deserialize;
use tracing::debug;

use crate::config::Config;

#[derive(Debug, Clone)]
pub struct FinditHit {
    pub bytes: Vec<u8>,
    pub source: String,
    pub doi: Option<String>,
}

#[derive(Debug, thiserror::Error)]
pub enum FinditError {
    #[error("findit disabled (FINDIT_URL unset)")]
    Disabled,
    #[error("findit returned no url ({0})")]
    NoUrl(String),
    #[error("{0}")]
    Other(String),
}

#[derive(Debug, Deserialize)]
struct FinditResp {
    #[allow(dead_code)]
    pmid: String,
    url: Option<String>,
    reason: Option<String>,
    doi: Option<String>,
    #[allow(dead_code)]
    journal: Option<String>,
    source: Option<String>,
}

pub async fn fetch_pdf(
    client: &Arc<Client>,
    cfg: &Config,
    pmid: &str,
) -> Result<FinditHit, FinditError> {
    let base = cfg.findit_url.as_deref().ok_or(FinditError::Disabled)?;
    let api = format!("{base}/findit/{pmid}");

    let resp = client
        .get(&api)
        .timeout(cfg.findit_timeout)
        .send()
        .await
        .map_err(|e| FinditError::Other(format!("findit request: {e}")))?;
    if !resp.status().is_success() {
        return Err(FinditError::Other(format!(
            "findit status {}",
            resp.status()
        )));
    }
    let body: FinditResp = resp
        .json()
        .await
        .map_err(|e| FinditError::Other(format!("findit decode: {e}")))?;

    let pdf_url = match body.url {
        Some(u) if !u.is_empty() => u,
        _ => {
            return Err(FinditError::NoUrl(
                body.reason.unwrap_or_else(|| "no url".into()),
            ));
        }
    };

    debug!(%pmid, %pdf_url, "findit resolved");

    let pdf_resp = client
        .get(&pdf_url)
        .timeout(cfg.download_timeout)
        .send()
        .await
        .map_err(|e| FinditError::Other(format!("pdf get: {e}")))?;
    if !pdf_resp.status().is_success() {
        return Err(FinditError::Other(format!(
            "pdf status {} ({pdf_url})",
            pdf_resp.status()
        )));
    }
    let ct = pdf_resp
        .headers()
        .get(reqwest::header::CONTENT_TYPE)
        .and_then(|v| v.to_str().ok())
        .unwrap_or("")
        .to_string();
    let bytes = pdf_resp
        .bytes()
        .await
        .map_err(|e| FinditError::Other(format!("pdf body: {e}")))?
        .to_vec();

    if !bytes.starts_with(b"%PDF-") && !ct.contains("application/pdf") {
        return Err(FinditError::Other(format!(
            "non-PDF payload (ct={ct}, {} bytes, url={pdf_url})",
            bytes.len()
        )));
    }

    let label = body.source.unwrap_or_else(|| "other".into());
    Ok(FinditHit {
        bytes,
        source: format!("findit:{label}"),
        doi: body.doi,
    })
}
