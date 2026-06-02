use std::sync::Arc;

use reqwest::{Client, Proxy};

use crate::config::Config;

const USER_AGENT: &str = "Mozilla/5.0 (Windows NT 10.0; Win64; x64) \
    AppleWebKit/537.36 (KHTML, like Gecko) Chrome/125.0.0.0 Safari/537.36";

pub fn build(cfg: &Config) -> anyhow::Result<Arc<Client>> {
    let mut builder = Client::builder()
        .user_agent(USER_AGENT)
        .timeout(cfg.http_timeout)
        .cookie_store(true)
        .danger_accept_invalid_certs(true)
        .gzip(true);

    if let Some(p) = &cfg.http_proxy {
        builder = builder.proxy(Proxy::http(p)?);
    }
    if let Some(p) = &cfg.https_proxy {
        builder = builder.proxy(Proxy::https(p)?);
    }

    Ok(Arc::new(builder.build()?))
}
