FROM ubuntu:24.04
ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y \
    curl \
    ca-certificates \
    bzip2 \
    git \
    && rm -rf /var/lib/apt/lists/*

# 1. Install the Pixi engine package manager (standalone binary).
RUN curl -fsSL --retry 5 --retry-all-errors --retry-delay 10 \
      -o /tmp/pixi.tar.gz \
      https://github.com/prefix-dev/pixi/releases/download/v0.58.0/pixi-x86_64-unknown-linux-musl.tar.gz \
    && tar -xzf /tmp/pixi.tar.gz -C /usr/local/bin pixi \
    && rm /tmp/pixi.tar.gz

# 2. Workspace destination
WORKDIR /comebin_env

# 3. Pixi manifest.
RUN cat <<'EOF' > pixi.toml
[workspace]
name = "comebin-gpu"
version = "1.0.0"
platforms = ["linux-64"]
channels = ["conda-forge", "bioconda"]

[system-requirements]
cuda = "12"

[dependencies]
python = "==3.11"
comebin = "*"
EOF

# 4. Solve and install the environment at BUILD time.
ENV CONDA_OVERRIDE_CUDA=12.0
RUN pixi install && rm -rf /root/.cache/rattler

# 5. Confirm the environment actually persisted into the image layer.
RUN ls -la /comebin_env/.pixi/envs/default/bin/python

# 6. Sanity check: fail the build if the CUDA torch build didn't land.
RUN pixi run python -c "import torch, sys; print('torch', torch.__version__, 'cuda', torch.version.cuda); sys.exit(1) if torch.version.cuda is None else None"

# 7. Put the environment on PATH so run_comebin.sh and friends resolve
#    without the `pixi run` wrapper.
ENV PATH="/comebin_env/.pixi/envs/default/bin:${PATH}"

ENTRYPOINT []
CMD ["/bin/bash"]
