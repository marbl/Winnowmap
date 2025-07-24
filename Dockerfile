FROM debian:stable AS builder

RUN apt-get update && apt-get install -y \
    git \
    build-essential \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-openssl-dev \
    ca-certificates \
    libomp-dev

COPY . /usr/src
WORKDIR /usr/src
RUN make -j8

FROM debian:stable-slim

RUN apt-get update && apt-get install -y \
    procps \
    libgomp1

COPY --from=builder /usr/src/bin/* /usr/local/bin/
