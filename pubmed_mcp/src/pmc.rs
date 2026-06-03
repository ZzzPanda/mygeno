use std::sync::Arc;

use reqwest::Client;
use serde_json::Value;
use tracing::debug;

use crate::config::Config;

#[derive(Debug, Clone)]
pub struct PmcHit {
    pub bytes: Vec<u8>,
    pub source: String,
}

#[derive(Debug, thiserror::Error)]
pub enum PmcError {
    #[error("no PMCID for this PMID")]
    NoPmcid,
    #[error("{0}")]
    Fetch(String),
}

const EPMC_SEARCH: &str = "https://www.ebi.ac.uk/europepmc/webservices/rest/search";

pub async fn lookup_pmcid(client: &Arc<Client>, pmid: &str) -> anyhow::Result<Option<String>> {
    let q = format!("EXT_ID:{pmid} AND SRC:MED");
    let resp = client
        .get(EPMC_SEARCH)
        .query(&[("query", q.as_str()), ("format", "json"), ("resultType", "lite")])
        .send()
        .await?
        .error_for_status()?;
    let json: Value = resp.json().await?;
    Ok(json
        .pointer("/resultList/result/0/pmcid")
        .and_then(Value::as_str)
        .map(|s| s.to_string()))
}

pub async fn fetch_pdf(
    client: &Arc<Client>,
    cfg: &Config,
    pmid: &str,
) -> Result<PmcHit, PmcError> {
    let pmcid = lookup_pmcid(client, pmid)
        .await
        .map_err(|e| PmcError::Fetch(format!("epmc lookup: {e}")))?
        .ok_or(PmcError::NoPmcid)?;

    let url = format!("https://europepmc.org/articles/{pmcid}?pdf=render");
    debug!(%pmid, %pmcid, %url, "trying Europe PMC");

    let resp = client
        .get(&url)
        .timeout(cfg.download_timeout)
        .send()
        .await
        .map_err(|e| PmcError::Fetch(format!("epmc get: {e}")))?;

    if !resp.status().is_success() {
        return Err(PmcError::Fetch(format!("epmc status {}", resp.status())));
    }
    let ct = resp
        .headers()
        .get(reqwest::header::CONTENT_TYPE)
        .and_then(|v| v.to_str().ok())
        .unwrap_or("")
        .to_string();
    let bytes = resp
        .bytes()
        .await
        .map_err(|e| PmcError::Fetch(format!("epmc body: {e}")))?
        .to_vec();

    if !bytes.starts_with(b"%PDF-") && !ct.contains("application/pdf") {
        return Err(PmcError::Fetch(format!(
            "epmc non-PDF (ct={ct}, {} bytes)",
            bytes.len()
        )));
    }

    Ok(PmcHit {
        bytes,
        source: format!("europepmc:{pmcid}"),
    })
}
