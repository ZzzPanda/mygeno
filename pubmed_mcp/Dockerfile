FROM rust:1-bookworm AS builder
WORKDIR /src

RUN apt-get update && apt-get install -y --no-install-recommends \
        pkg-config libssl-dev \
    && rm -rf /var/lib/apt/lists/*

COPY Cargo.toml Cargo.lock* ./
RUN mkdir src && echo 'fn main() {}' > src/main.rs && \
    cargo build --release --bin pmid_downloader && \
    rm -rf src target/release/deps/pmid_downloader*

COPY migrations ./migrations
COPY src ./src
RUN cargo build --release --bin pmid_downloader

FROM debian:bookworm-slim
RUN apt-get update && apt-get install -y --no-install-recommends \
        ca-certificates \
    && rm -rf /var/lib/apt/lists/*

RUN useradd -u 1000 -m -s /bin/bash app && mkdir -p /data && chown app:app /data
USER app
WORKDIR /home/app

COPY --from=builder /src/target/release/pmid_downloader /usr/local/bin/pmid_downloader
COPY --from=builder /src/migrations /migrations

ENV SQLITE_PATH=/data/meta.db \
    MCP_BIND=0.0.0.0:8080 \
    MCP_PATH=/mcp \
    RUST_LOG=info,pmid_downloader=debug

EXPOSE 8080
ENTRYPOINT ["/usr/local/bin/pmid_downloader"]
