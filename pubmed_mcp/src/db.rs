use std::path::Path;

use sqlx::{
    sqlite::{SqliteConnectOptions, SqlitePoolOptions},
    SqlitePool,
};

#[derive(Debug, Clone, sqlx::FromRow)]
pub struct DownloadRow {
    pub pmid: String,
    pub doi: String,
    pub object_key: String,
    pub bucket: String,
    pub size_bytes: i64,
    pub sha256: String,
    pub source_mirror: String,
    pub downloaded_at: i64,
}

#[derive(Clone)]
pub struct Db {
    pool: SqlitePool,
}

impl Db {
    pub async fn open(path: &str) -> anyhow::Result<Self> {
        if let Some(parent) = Path::new(path).parent() {
            if !parent.as_os_str().is_empty() {
                std::fs::create_dir_all(parent).ok();
            }
        }
        let opts = SqliteConnectOptions::new()
            .filename(path)
            .create_if_missing(true)
            .foreign_keys(true);
        let pool = SqlitePoolOptions::new()
            .max_connections(5)
            .connect_with(opts)
            .await?;
        sqlx::migrate!("./migrations").run(&pool).await?;
        Ok(Self { pool })
    }

    pub async fn lookup(&self, pmid: &str) -> anyhow::Result<Option<DownloadRow>> {
        let row = sqlx::query_as::<_, DownloadRow>(
            "SELECT pmid, doi, object_key, bucket, size_bytes, sha256, source_mirror, downloaded_at \
             FROM downloads WHERE pmid = ?",
        )
        .bind(pmid)
        .fetch_optional(&self.pool)
        .await?;
        Ok(row)
    }

    pub async fn insert(&self, row: &DownloadRow) -> anyhow::Result<()> {
        sqlx::query(
            "INSERT OR REPLACE INTO downloads \
             (pmid, doi, object_key, bucket, size_bytes, sha256, source_mirror, downloaded_at) \
             VALUES (?, ?, ?, ?, ?, ?, ?, ?)",
        )
        .bind(&row.pmid)
        .bind(&row.doi)
        .bind(&row.object_key)
        .bind(&row.bucket)
        .bind(row.size_bytes)
        .bind(&row.sha256)
        .bind(&row.source_mirror)
        .bind(row.downloaded_at)
        .execute(&self.pool)
        .await?;
        Ok(())
    }

    pub async fn list(&self, limit: i64, offset: i64) -> anyhow::Result<Vec<DownloadRow>> {
        let rows = sqlx::query_as::<_, DownloadRow>(
            "SELECT pmid, doi, object_key, bucket, size_bytes, sha256, source_mirror, downloaded_at \
             FROM downloads ORDER BY downloaded_at DESC LIMIT ? OFFSET ?",
        )
        .bind(limit)
        .bind(offset)
        .fetch_all(&self.pool)
        .await?;
        Ok(rows)
    }
}
