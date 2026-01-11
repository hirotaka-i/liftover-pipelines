# Minimal Docker image for GWAS summary statistics liftover
# Includes bcftools with liftover plugin and Python dependencies

FROM ubuntu:22.04

ENV DEBIAN_FRONTEND=noninteractive
ENV PATH=/usr/local/bin:$PATH

# Install system dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    ca-certificates wget curl unzip \
    build-essential cmake gcc g++ make autoconf automake libtool \
    python3.10 python3-pip \
    zlib1g-dev libbz2-dev libgsl-dev \
    liblzma-dev libcurl4-openssl-dev \
    libssl-dev libncurses5-dev tabix \
    && rm -rf /var/lib/apt/lists/*

# Install Python packages
RUN pip3 install --no-cache-dir \
    pandas==2.1.4 \
    numpy==1.26.3

# ============================================
# HTSlib, BCFtools with liftover plugin
# ============================================

ENV HTSLIB_VERSION=1.20
ENV BCFTOOLS_VERSION=1.20

# Install HTSlib
RUN wget -q https://github.com/samtools/htslib/releases/download/${HTSLIB_VERSION}/htslib-${HTSLIB_VERSION}.tar.bz2 -O /tmp/htslib.tar.bz2 && \
    cd /tmp && tar xjf htslib.tar.bz2 && \
    cd htslib-${HTSLIB_VERSION} && \
    autoreconf -i && \
    ./configure && \
    make && make install && \
    ldconfig && \
    rm -rf /tmp/htslib*

# Install BCFtools with liftover plugin
RUN wget -q https://github.com/samtools/bcftools/releases/download/${BCFTOOLS_VERSION}/bcftools-${BCFTOOLS_VERSION}.tar.bz2 -O /tmp/bcftools.tar.bz2 && \
    cd /tmp && tar xjf bcftools.tar.bz2 && \
    rm bcftools.tar.bz2

# Download liftover plugin
RUN wget -q https://raw.githubusercontent.com/freeseek/score/a9ffa435913101439974fce7c4812235b61e6df5/liftover.c \
    -O /tmp/bcftools-${BCFTOOLS_VERSION}/plugins/liftover.c

# Compile and install BCFtools
RUN cd /tmp/bcftools-${BCFTOOLS_VERSION} && \
    ./configure --enable-libgsl --with-htslib=system && \
    make && make install && \
    rm -rf /tmp/bcftools-${BCFTOOLS_VERSION}

# ============================================
# UCSC liftOver tool
# ============================================

RUN wget -q http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/liftOver -O /usr/local/bin/liftOver && \
    chmod +x /usr/local/bin/liftOver

# Create directories for reference files
RUN mkdir -p /workspace /references

# Copy liftover script (simple version only)
COPY codes/liftover_sumstats_simple.py /usr/local/bin/
COPY liftover /usr/local/bin/
RUN chmod +x /usr/local/bin/liftover_sumstats_simple.py /usr/local/bin/liftover

WORKDIR /workspace

# Set the entrypoint to our wrapper script
ENTRYPOINT ["/usr/local/bin/liftover"]

# Verify installation
RUN python3 --version && bcftools --version && liftOver 2>&1 | head -1
