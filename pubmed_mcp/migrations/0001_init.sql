CREATE TABLE IF NOT EXISTS downloads (
    pmid          TEXT PRIMARY KEY,
    doi           TEXT NOT NULL,
    object_key    TEXT NOT NULL,
    bucket        TEXT NOT NULL,
    size_bytes    INTEGER NOT NULL,
    sha256        TEXT NOT NULL,
    source_mirror TEXT NOT NULL,
    downloaded_at INTEGER NOT NULL
);

CREATE INDEX IF NOT EXISTS idx_downloads_doi ON downloads(doi);
CREATE INDEX IF NOT EXISTS idx_downloads_at  ON downloads(downloaded_at);
