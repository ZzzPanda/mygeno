use rmcp::{
    handler::server::{router::tool::ToolRouter, wrapper::Parameters},
    model::{CallToolResult, Content, Implementation, ProtocolVersion, ServerCapabilities, ServerInfo},
    schemars, tool, tool_handler, tool_router, ErrorData as McpError, ServerHandler,
};
use serde::{Deserialize, Serialize};
use serde_json::json;

use crate::pipeline::Pipeline;

#[derive(Debug, Deserialize, schemars::JsonSchema)]
pub struct DownloadPmidArgs {
    /// PubMed ID (digits only)
    pub pmid: String,
}

#[derive(Debug, Deserialize, schemars::JsonSchema)]
pub struct DownloadPmidsArgs {
    /// List of PubMed IDs to download in sequence
    pub pmids: Vec<String>,
}

#[derive(Debug, Deserialize, schemars::JsonSchema)]
pub struct GetPdfUriArgs {
    /// PubMed ID to look up
    pub pmid: String,
}

#[derive(Debug, Deserialize, schemars::JsonSchema)]
pub struct ListDownloadedArgs {
    /// Max records to return (default 50)
    #[serde(default)]
    pub limit: Option<i64>,
    /// Records to skip (default 0)
    #[serde(default)]
    pub offset: Option<i64>,
}

#[derive(Debug, Serialize, schemars::JsonSchema)]
pub struct PdfRecord {
    pub pmid: String,
    pub doi: String,
    pub uri: String,
    pub size_bytes: i64,
    pub downloaded_at: i64,
}

#[derive(Clone)]
pub struct PmidServer {
    pipeline: Pipeline,
    tool_router: ToolRouter<PmidServer>,
}

#[tool_router]
impl PmidServer {
    pub fn new(pipeline: Pipeline) -> Self {
        Self {
            pipeline,
            tool_router: Self::tool_router(),
        }
    }

    #[tool(description = "Download PDF for a single PubMed ID. Returns S3 URI and metadata; cached on repeat.")]
    async fn download_pmid(
        &self,
        Parameters(args): Parameters<DownloadPmidArgs>,
    ) -> Result<CallToolResult, McpError> {
        let outcome = self.pipeline.download_one(&args.pmid).await;
        Ok(CallToolResult::success(vec![Content::text(
            serde_json::to_string(&outcome).unwrap_or_else(|_| "{}".into()),
        )]))
    }

    #[tool(description = "Download PDFs for multiple PubMed IDs sequentially with human-paced delays.")]
    async fn download_pmids(
        &self,
        Parameters(args): Parameters<DownloadPmidsArgs>,
    ) -> Result<CallToolResult, McpError> {
        let results = self.pipeline.download_many(&args.pmids).await;
        Ok(CallToolResult::success(vec![Content::text(
            serde_json::to_string(&results).unwrap_or_else(|_| "[]".into()),
        )]))
    }

    #[tool(description = "Look up the S3 URI of an already-downloaded PMID without re-fetching.")]
    async fn get_pdf_uri(
        &self,
        Parameters(args): Parameters<GetPdfUriArgs>,
    ) -> Result<CallToolResult, McpError> {
        let pmid = args.pmid.trim().to_string();
        let row = self
            .pipeline
            .db
            .lookup(&pmid)
            .await
            .map_err(|e| McpError::internal_error(e.to_string(), None))?;

        let payload = match row {
            Some(r) => json!({
                "pmid": r.pmid,
                "uri": self.pipeline.storage.uri_for(&r.object_key),
                "status": "found",
            }),
            None => json!({"pmid": pmid, "status": "missing"}),
        };
        Ok(CallToolResult::success(vec![Content::text(
            payload.to_string(),
        )]))
    }

    #[tool(description = "List downloaded PDF records ordered by most recent first.")]
    async fn list_downloaded(
        &self,
        Parameters(args): Parameters<ListDownloadedArgs>,
    ) -> Result<CallToolResult, McpError> {
        let limit = args.limit.unwrap_or(50).clamp(1, 500);
        let offset = args.offset.unwrap_or(0).max(0);
        let rows = self
            .pipeline
            .db
            .list(limit, offset)
            .await
            .map_err(|e| McpError::internal_error(e.to_string(), None))?;
        let records: Vec<PdfRecord> = rows
            .into_iter()
            .map(|r| PdfRecord {
                uri: self.pipeline.storage.uri_for(&r.object_key),
                pmid: r.pmid,
                doi: r.doi,
                size_bytes: r.size_bytes,
                downloaded_at: r.downloaded_at,
            })
            .collect();
        Ok(CallToolResult::success(vec![Content::text(
            serde_json::to_string(&records).unwrap_or_else(|_| "[]".into()),
        )]))
    }
}

#[tool_handler]
impl ServerHandler for PmidServer {
    fn get_info(&self) -> ServerInfo {
        ServerInfo::new(
            ServerCapabilities::builder().enable_tools().build(),
        )
        .with_server_info(Implementation::from_build_env())
        .with_protocol_version(ProtocolVersion::V_2024_11_05)
        .with_instructions(
            "PMID PDF downloader. Tools: download_pmid, download_pmids, get_pdf_uri, list_downloaded. \
             PDFs land in MinIO (S3-compatible); metadata in SQLite. Tools return s3:// URIs."
                .to_string(),
        )
    }
}
