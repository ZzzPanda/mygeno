use std::sync::Arc;

use reqwest::Client;
use serde_json::Value;

const ESUMMARY_URL: &str = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi";

pub async fn fetch_doi(client: &Arc<Client>, pmid: &str) -> anyhow::Result<Option<String>> {
    let resp = client
        .get(ESUMMARY_URL)
        .query(&[("db", "pubmed"), ("id", pmid), ("retmode", "json")])
        .send()
        .await?
        .error_for_status()?;

    let json: Value = resp.json().await?;
    let article = json
        .get("result")
        .and_then(|r| r.get(pmid))
        .cloned()
        .unwrap_or(Value::Null);

    let mut doi = article
        .get("elocationid")
        .and_then(Value::as_str)
        .unwrap_or("")
        .to_string();
    if doi.is_empty() {
        doi = article
            .get("doi")
            .and_then(Value::as_str)
            .unwrap_or("")
            .to_string();
    }

    let doi = doi.trim();
    let doi = doi.strip_prefix("doi:").unwrap_or(doi).trim().to_string();
    if doi.is_empty() {
        return Ok(None);
    }
    Ok(Some(doi))
}
