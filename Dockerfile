FROM debian:bookworm-slim

RUN apt-get update && \
    apt-get install -y --no-install-recommends procps ca-certificates && \
    update-ca-certificates && \
    rm -rf /var/lib/apt /var/lib/dpkg /var/lib/cache /var/lib/log

WORKDIR /app

# The release workflow supplies statically linked musl binaries so the image does not depend on
# the build runner's glibc version.
COPY target/x86_64-unknown-linux-musl/release/sage /app/sage

COPY THIRD_PARTY_NOTICES.md /app/licenses/THIRD_PARTY_NOTICES.md
COPY vendor/filemanager/LICENSE /app/licenses/filemanager.txt
COPY vendor/filemanager/LICENSE-APACHE /app/licenses/filemanager-APACHE.txt

ENV PATH="/app:$PATH"
