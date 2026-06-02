use std::sync::Arc;

use aws_credential_types::Credentials;
use aws_sdk_s3::{
    config::{BehaviorVersion, Region},
    primitives::ByteStream,
    Client,
};

use crate::config::Config;

#[derive(Clone)]
pub struct Storage {
    client: Arc<Client>,
    bucket: String,
}

impl Storage {
    pub async fn from_config(cfg: &Config) -> anyhow::Result<Self> {
        let creds = Credentials::new(
            &cfg.s3_access_key,
            &cfg.s3_secret_key,
            None,
            None,
            "static",
        );
        let s3_cfg = aws_sdk_s3::Config::builder()
            .behavior_version(BehaviorVersion::latest())
            .region(Region::new(cfg.s3_region.clone()))
            .endpoint_url(&cfg.s3_endpoint)
            .credentials_provider(creds)
            .force_path_style(cfg.s3_force_path_style)
            .build();
        let client = Client::from_conf(s3_cfg);
        Ok(Self {
            client: Arc::new(client),
            bucket: cfg.s3_bucket.clone(),
        })
    }

    pub fn bucket(&self) -> &str {
        &self.bucket
    }

    pub fn uri_for(&self, key: &str) -> String {
        format!("s3://{}/{}", self.bucket, key)
    }

    pub async fn put(&self, key: &str, bytes: Vec<u8>) -> anyhow::Result<()> {
        self.client
            .put_object()
            .bucket(&self.bucket)
            .key(key)
            .content_type("application/pdf")
            .body(ByteStream::from(bytes))
            .send()
            .await?;
        Ok(())
    }

    pub async fn ensure_bucket(&self) -> anyhow::Result<()> {
        match self.client.head_bucket().bucket(&self.bucket).send().await {
            Ok(_) => Ok(()),
            Err(_) => {
                self.client
                    .create_bucket()
                    .bucket(&self.bucket)
                    .send()
                    .await?;
                Ok(())
            }
        }
    }
}
