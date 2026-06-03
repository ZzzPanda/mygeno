use std::sync::Arc;

use rmcp::transport::streamable_http_server::{
    session::local::LocalSessionManager, StreamableHttpServerConfig, StreamableHttpService,
};
use tokio_util::sync::CancellationToken;
use tracing::info;
use tracing_subscriber::{layer::SubscriberExt, util::SubscriberInitExt};

mod config;
mod db;
mod findit;
mod http_client;
mod pipeline;
mod pmc;
mod pubmed;
mod scihub;
mod server;
mod storage;

use crate::{
    config::Config, db::Db, pipeline::Pipeline, server::PmidServer, storage::Storage,
};

#[tokio::main]
async fn main() -> anyhow::Result<()> {
    tracing_subscriber::registry()
        .with(
            tracing_subscriber::EnvFilter::try_from_default_env()
                .unwrap_or_else(|_| "info,pmid_downloader=debug".into()),
        )
        .with(tracing_subscriber::fmt::layer())
        .init();

    let cfg = Config::from_env()?;
    info!(?cfg.bind, ?cfg.mcp_path, ?cfg.s3_endpoint, ?cfg.s3_bucket, ?cfg.sqlite_path, "starting pmid_downloader");

    let http = http_client::build(&cfg)?;
    let db = Db::open(&cfg.sqlite_path).await?;
    let storage = Storage::from_config(&cfg).await?;
    if let Err(e) = storage.ensure_bucket().await {
        tracing::warn!(error = %e, "ensure_bucket failed (may already exist via init container)");
    }

    scihub::warmup(&http, &cfg).await;

    let pipeline = Pipeline {
        cfg: cfg.clone(),
        http,
        db,
        storage,
    };
    let pipeline = Arc::new(pipeline);

    let ct = CancellationToken::new();
    let service = StreamableHttpService::new(
        {
            let pipeline = pipeline.clone();
            move || Ok(PmidServer::new((*pipeline).clone()))
        },
        LocalSessionManager::default().into(),
        StreamableHttpServerConfig::default().with_cancellation_token(ct.child_token()),
    );

    let router = axum::Router::new().nest_service(&cfg.mcp_path, service);
    let listener = tokio::net::TcpListener::bind(&cfg.bind).await?;
    info!(bind = %cfg.bind, path = %cfg.mcp_path, "MCP listening");

    axum::serve(listener, router)
        .with_graceful_shutdown(async move {
            tokio::signal::ctrl_c().await.ok();
            ct.cancel();
        })
        .await?;
    Ok(())
}
