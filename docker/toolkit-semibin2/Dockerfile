FROM ubuntu:24.04 AS builder
ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y \
    curl \
    ca-certificates \
    bzip2 \
    git \
    && rm -rf /var/lib/apt/lists/*

RUN curl -fsSL --retry 5 --retry-all-errors --retry-delay 10 \
      -o /tmp/pixi.tar.gz \
      https://github.com/prefix-dev/pixi/releases/download/v0.58.0/pixi-x86_64-unknown-linux-musl.tar.gz \
    && tar -xzf /tmp/pixi.tar.gz -C /usr/local/bin pixi \
    && rm /tmp/pixi.tar.gz

WORKDIR /semibin_env

RUN cat <<'EOF' > pixi.toml
[workspace]
name = "semibin-gpu"
version = "1.0.0"
platforms = ["linux-64"]
channels = ["conda-forge", "bioconda"]

[system-requirements]
cuda = "12.0"

[dependencies]
semibin = ">=2.4.1,<3"
pytorch-gpu = "*"
bedtools = "*"
hmmer = "*"
samtools = "*"
prodigal = "*"
EOF

# No GPU in the builder, so mock the __cuda virtual package.
ENV CONDA_OVERRIDE_CUDA=12.0
RUN pixi install

# Trim runtime-unnecessary files. Headers and static libs are a large part
# of the conda CUDA stack and nothing here compiles at runtime.
RUN find /semibin_env/.pixi/envs/default -type d -name include -prune -exec rm -rf {} + 2>/dev/null || true \
    && find /semibin_env/.pixi/envs/default -name '*.a' -delete 2>/dev/null || true \
    && find /semibin_env/.pixi/envs/default -type d -name '__pycache__' -prune -exec rm -rf {} + 2>/dev/null || true

FROM ubuntu:24.04
ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y \
    ca-certificates \
    bzip2 \
    && rm -rf /var/lib/apt/lists/*

COPY --from=builder /semibin_env /semibin_env

ENV PATH="/semibin_env/.pixi/envs/default/bin:${PATH}"

RUN python -c "import torch, sys; print('torch', torch.__version__, 'cuda', torch.version.cuda); sys.exit(1) if torch.version.cuda is None else None"
RUN SemiBin2 --version

ENTRYPOINT []
CMD ["/bin/bash"]
