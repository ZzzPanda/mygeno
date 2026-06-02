use std::{sync::Arc, time::Duration};

use rand::Rng;
use regex::Regex;
use reqwest::{Client, Url};
use scraper::{Html, Selector};
use tracing::{debug, warn};

use crate::config::Config;

#[derive(Debug, Clone)]
pub struct PdfHit {
    pub bytes: Vec<u8>,
    pub mirror: String,
    pub source_url: String,
}

#[derive(Debug, thiserror::Error)]
pub enum ScihubError {
    #[error("all mirrors exhausted: {0}")]
    Exhausted(String),
}

pub async fn warmup(client: &Arc<Client>, cfg: &Config) {
    if let Some(mirror) = cfg.mirrors.first() {
        let url = format!("http://{mirror}/");
        if let Err(e) = client.get(&url).send().await {
            warn!(error = %e, %url, "warmup failed");
        }
        sleep(cfg.human_delay, 2.0..5.0).await;
    }
}

pub async fn fetch_pdf(
    client: &Arc<Client>,
    cfg: &Config,
    doi: &str,
) -> Result<PdfHit, ScihubError> {
    let mut errs = Vec::new();
    for (i, mirror) in cfg.mirrors.iter().enumerate() {
        if i > 0 {
            sleep(cfg.human_delay, 3.0..8.0).await;
        }
        match try_mirror(client, cfg, mirror, doi).await {
            Ok(hit) => return Ok(hit),
            Err(e) => {
                debug!(%mirror, error = %e, "mirror failed");
                errs.push(format!("{mirror}: {e}"));
            }
        }
    }
    Err(ScihubError::Exhausted(errs.join("; ")))
}

async fn try_mirror(
    client: &Arc<Client>,
    cfg: &Config,
    mirror: &str,
    doi: &str,
) -> anyhow::Result<PdfHit> {
    for (i, proto) in ["http", "https"].iter().enumerate() {
        if i > 0 {
            sleep(cfg.human_delay, 1.0..2.5).await;
        }
        let scihub_url = format!("{proto}://{mirror}/{doi}");

        let resp = match client.get(&scihub_url).send().await {
            Ok(r) => r,
            Err(e) => {
                debug!(%scihub_url, error = %e, "request failed");
                continue;
            }
        };
        if !resp.status().is_success() {
            debug!(%scihub_url, status = %resp.status(), "non-200");
            continue;
        }

        sleep(cfg.human_delay, 2.0..5.0).await;

        let final_url = resp.url().clone();
        let content_type = resp
            .headers()
            .get(reqwest::header::CONTENT_TYPE)
            .and_then(|v| v.to_str().ok())
            .unwrap_or("")
            .to_string();
        let bytes = resp.bytes().await?;

        if bytes.starts_with(b"%PDF-") || content_type.contains("application/pdf") {
            return Ok(PdfHit {
                bytes: bytes.to_vec(),
                mirror: mirror.to_string(),
                source_url: scihub_url,
            });
        }

        let html = String::from_utf8_lossy(&bytes).to_string();
        if is_captcha_or_blocked(&html) {
            debug!(%scihub_url, len = html.len(), "captcha/blocked");
            continue;
        }

        let pdf_url = match find_pdf_url(&html, &final_url) {
            Some(u) => u,
            None => {
                debug!(%scihub_url, len = html.len(), "no pdf url found");
                continue;
            }
        };

        sleep(cfg.human_delay, 0.8..2.5).await;

        let pdf_resp = client
            .get(pdf_url.clone())
            .timeout(cfg.download_timeout)
            .header(reqwest::header::REFERER, &scihub_url)
            .send()
            .await?;
        if !pdf_resp.status().is_success() {
            debug!(%pdf_url, status = %pdf_resp.status(), "pdf non-200");
            continue;
        }
        let pdf_ct = pdf_resp
            .headers()
            .get(reqwest::header::CONTENT_TYPE)
            .and_then(|v| v.to_str().ok())
            .unwrap_or("")
            .to_string();
        let pdf_bytes = pdf_resp.bytes().await?;
        if pdf_bytes.starts_with(b"%PDF-") || pdf_ct.contains("application/pdf") {
            return Ok(PdfHit {
                bytes: pdf_bytes.to_vec(),
                mirror: mirror.to_string(),
                source_url: pdf_url,
            });
        }
        debug!(%pdf_url, len = pdf_bytes.len(), "non-pdf payload");
    }
    anyhow::bail!("all protocols failed");
}

pub fn validate_pdf(bytes: &[u8]) -> Result<(), &'static str> {
    if bytes.len() < 1000 {
        return Err("file too small");
    }
    if !bytes.starts_with(b"%PDF-") {
        return Err("not a PDF");
    }
    let tail_start = bytes.len().saturating_sub(1024);
    if !bytes[tail_start..].windows(5).any(|w| w == b"%%EOF") {
        return Err("PDF incomplete (no %%EOF)");
    }
    Ok(())
}

fn is_captcha_or_blocked(html: &str) -> bool {
    let body = html.trim();
    if body.len() < 200 {
        return true;
    }
    if body.contains("/storage/")
        || body.contains("citation_pdf_url")
        || body.contains("pdf2md")
    {
        return false;
    }
    let doc = Html::parse_document(body);
    if let Ok(sel) = Selector::parse(r#"object[type="application/pdf"]"#) {
        if doc.select(&sel).next().is_some() {
            return false;
        }
    }
    if let Ok(sel) = Selector::parse("#cf-wrapper, .cf-browser-verification") {
        if doc.select(&sel).next().is_some() {
            return true;
        }
    }
    let lower = body.to_lowercase();
    for kw in ["recaptcha", "cf-captcha", "ddos protection"] {
        if lower.contains(kw) {
            return true;
        }
    }
    false
}

fn find_pdf_url(html: &str, page_url: &Url) -> Option<String> {
    let doc = Html::parse_document(html);

    if let Ok(sel) = Selector::parse(r#"meta[name="citation_pdf_url"]"#) {
        if let Some(el) = doc.select(&sel).next() {
            if let Some(c) = el.value().attr("content") {
                return Some(absolutize(c, page_url));
            }
        }
    }
    if let Ok(sel) = Selector::parse(r#"object[type="application/pdf"]"#) {
        if let Some(el) = doc.select(&sel).next() {
            if let Some(d) = el.value().attr("data") {
                return Some(absolutize(d, page_url));
            }
        }
    }
    if let Ok(sel) = Selector::parse("div.download a[href]") {
        if let Some(el) = doc.select(&sel).next() {
            if let Some(h) = el.value().attr("href") {
                return Some(absolutize(h, page_url));
            }
        }
    }
    if let Ok(sel) = Selector::parse("iframe") {
        for el in doc.select(&sel) {
            if let Some(s) = el.value().attr("src") {
                if !s.starts_with("about:") {
                    return Some(absolutize(s, page_url));
                }
            }
        }
    }
    if let Ok(sel) = Selector::parse("embed") {
        if let Some(el) = doc.select(&sel).next() {
            if let Some(s) = el.value().attr("src") {
                return Some(absolutize(s, page_url));
            }
        }
    }
    if let Ok(sel) = Selector::parse("button") {
        let onclick_re =
            Regex::new(r#"location\s*\.?\s*href\s*=\s*['"]([^'"]+)['"]"#).ok()?;
        for el in doc.select(&sel) {
            if let Some(oc) = el.value().attr("onclick") {
                if let Some(m) = onclick_re.captures(oc).and_then(|c| c.get(1)) {
                    return Some(absolutize(m.as_str(), page_url));
                }
            }
        }
    }
    if let Ok(re) = Regex::new(r#"(/storage/[^\s"']+\.pdf)"#) {
        if let Some(m) = re.captures(html).and_then(|c| c.get(1)) {
            return Some(absolutize(m.as_str(), page_url));
        }
    }
    if let Ok(sel) = Selector::parse("a[href]") {
        for el in doc.select(&sel) {
            if let Some(h) = el.value().attr("href") {
                if h.to_lowercase().contains(".pdf") {
                    return Some(absolutize(h, page_url));
                }
            }
        }
    }
    None
}

fn absolutize(href: &str, base: &Url) -> String {
    if let Some(rest) = href.strip_prefix("//") {
        return format!("https://{rest}");
    }
    if href.starts_with("http://") || href.starts_with("https://") {
        return href.to_string();
    }
    base.join(href).map(|u| u.to_string()).unwrap_or_else(|_| href.to_string())
}

async fn sleep(enabled: bool, range: std::ops::Range<f64>) {
    if !enabled {
        return;
    }
    let secs = rand::thread_rng().gen_range(range);
    tokio::time::sleep(Duration::from_secs_f64(secs)).await;
}
