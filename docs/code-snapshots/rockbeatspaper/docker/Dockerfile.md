---
source: /Users/byron/projects/rockbeatspaper/docker/Dockerfile
relative: rockbeatspaper/docker/Dockerfile
generated_at: 2025-12-23 10:28
---

```docker
FROM debian:12

WORKDIR /app

RUN DEBIAN_FRONTEND=noninteractive apt-get update && \
    DEBIAN_FRONTEND=noninteractive apt-get install --yes curl wget && \
    curl -fsSL https://deb.nodesource.com/setup_20.x | bash && \
    DEBIAN_FRONTEND=noninteractive apt-get install --yes nodejs && \
    DEBIAN_FRONTEND=noninteractive apt-get install --yes wget build-essential libcairo2-dev libpango1.0-dev libjpeg-dev libgif-dev librsvg2-dev vim less iputils-ping sudo libsecret-1-0 command-not-found rsync man-db netcat-openbsd dnsutils procps lsof tini && \
    DEBIAN_FRONTEND=noninteractive apt-get update

# Install Python 3.12 from source (required by thrml-hack package)
RUN DEBIAN_FRONTEND=noninteractive apt-get install --yes \
    libssl-dev zlib1g-dev libncurses5-dev libncursesw5-dev \
    libreadline-dev libsqlite3-dev libgdbm-dev libdb5.3-dev \
    libbz2-dev libexpat1-dev liblzma-dev libffi-dev uuid-dev && \
    cd /tmp && \
    wget https://www.python.org/ftp/python/3.12.8/Python-3.12.8.tgz && \
    tar xzf Python-3.12.8.tgz && \
    cd Python-3.12.8 && \
    ./configure --enable-optimizations --with-ensurepip=install && \
    make -j$(nproc) && \
    make altinstall && \
    ln -sf /usr/local/bin/python3.12 /usr/bin/python3 && \
    ln -sf /usr/local/bin/pip3.12 /usr/bin/pip3 && \
    cd / && \
    rm -rf /tmp/Python-3.12.8 /tmp/Python-3.12.8.tgz

RUN npm install -g pnpm
RUN pnpm set store-dir /app/node_modules/.pnpm-store

```
