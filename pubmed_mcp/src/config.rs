use std::{env, time::Duration};

#[derive(Debug, Clone)]
pub struct Config {
    pub bind: String,
    pub mcp_path: String,

    pub s3_endpoint: String,
    pub s3_region: String,
    pub s3_bucket: String,
    pub s3_access_key: String,
    pub s3_secret_key: String,
    pub s3_force_path_style: bool,

    pub sqlite_path: String,

    pub mirrors: Vec<String>,
    pub http_timeout: Duration,
    pub download_timeout: Duration,
    pub human_delay: bool,

    pub findit_url: Option<String>,
    pub findit_timeout: Duration,

    pub http_proxy: Option<String>,
    pub https_proxy: Option<String>,
}

impl Config {
    pub fn from_env() -> anyhow::Result<Self> {
        Ok(Self {
            bind: env::var("MCP_BIND").unwrap_or_else(|_| "0.0.0.0:8080".into()),
            mcp_path: env::var("MCP_PATH").unwrap_or_else(|_| "/mcp".into()),

            s3_endpoint: env::var("S3_ENDPOINT")
                .unwrap_or_else(|_| "http://minio:9000".into()),
            s3_region: env::var("S3_REGION").unwrap_or_else(|_| "us-east-1".into()),
            s3_bucket: env::var("S3_BUCKET").unwrap_or_else(|_| "pmid-pdfs".into()),
            s3_access_key: env::var("S3_ACCESS_KEY").unwrap_or_else(|_| "minio".into()),
            s3_secret_key: env::var("S3_SECRET_KEY").unwrap_or_else(|_| "minio12345".into()),
            s3_force_path_style: env::var("S3_FORCE_PATH_STYLE")
                .map(|v| v == "true" || v == "1")
                .unwrap_or(true),

            sqlite_path: env::var("SQLITE_PATH").unwrap_or_else(|_| "/data/meta.db".into()),

            mirrors: env::var("SCIHUB_MIRRORS")
                .unwrap_or_else(|_| "sci-hub.ru,sci-hub.st,sci-hub.se".into())
                .split(',')
                .map(|s| s.trim().to_string())
                .filter(|s| !s.is_empty())
                .collect(),
            http_timeout: Duration::from_secs(
                env::var("HTTP_TIMEOUT_SECS")
                    .ok()
                    .and_then(|s| s.parse().ok())
                    .unwrap_or(30),
            ),
            download_timeout: Duration::from_secs(
                env::var("DOWNLOAD_TIMEOUT_SECS")
                    .ok()
                    .and_then(|s| s.parse().ok())
                    .unwrap_or(120),
            ),
            human_delay: env::var("HUMAN_DELAY")
                .map(|v| v != "false" && v != "0")
                .unwrap_or(true),

            findit_url: env::var("FINDIT_URL")
                .ok()
                .map(|s| s.trim_end_matches('/').to_string())
                .filter(|s| !s.is_empty()),
            findit_timeout: Duration::from_secs(
                env::var("FINDIT_TIMEOUT_SECS")
                    .ok()
                    .and_then(|s| s.parse().ok())
                    .unwrap_or(30),
            ),

            http_proxy: env::var("HTTP_PROXY").ok().filter(|s| !s.is_empty()),
            https_proxy: env::var("HTTPS_PROXY").ok().filter(|s| !s.is_empty()),
        })
    }
}
