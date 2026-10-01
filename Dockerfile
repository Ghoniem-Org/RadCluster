# RadCluster — reproducible image for the EUROFER-97 cluster-dynamics suite.
#
# Two stages.  The first compiles the C++ solver against SUNDIALS/CVODE; the
# second carries only the runtime, so the published image does not ship a
# compiler toolchain or the SUNDIALS headers.
#
# Build:
#   docker build -t radcluster .
# Run a shell with the module on PYTHONPATH:
#   docker run --rm -it radcluster
# Run the solver directly:
#   docker run --rm radcluster solver --help
# Mount a host directory to keep the timestamped run output:
#   docker run --rm -v "$PWD/output:/work/output" radcluster \
#       python -c "from py_utils.simulation import ..."

# ─────────────────────────────────────────────────────────────────────────────
# Stage 1 — build the solver
# ─────────────────────────────────────────────────────────────────────────────
FROM ubuntu:24.04 AS build

ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get install -y --no-install-recommends \
        build-essential \
        cmake \
        libsundials-dev \
        liblapack-dev \
        libopenblas-dev \
        libsuitesparse-dev \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /src
COPY radcluster_code/cpp_utils ./cpp_utils

# -march=native would bake this builder's instruction set into a binary that
# other people pull and run; x86-64-v2 (SSE4.2 / POPCNT, 2009 and later) is the
# portable floor.  Fast-math stays off here for the same reason it is off in the
# local build: it perturbs IEEE-754 semantics and can make CVODE's error
# estimator diverge on these stiff systems.
RUN cmake -S cpp_utils -B build \
        -DCMAKE_BUILD_TYPE=Release \
        -DRADCLUSTER_ARCH_FLAGS=-march=x86-64-v2 \
    && cmake --build build --parallel "$(nproc)"

# The solver lands in build/ or build/Release/ depending on the generator.
RUN set -eux; \
    found="$(find build -name solver -type f -perm -u+x | head -1)"; \
    test -n "$found"; \
    install -Dm755 "$found" /out/solver; \
    /out/solver --help >/dev/null 2>&1 || true

# ─────────────────────────────────────────────────────────────────────────────
# Stage 2 — runtime
# ─────────────────────────────────────────────────────────────────────────────
FROM ubuntu:24.04

LABEL org.opencontainers.image.title="RadCluster" \
      org.opencontainers.image.description="Graph-based cluster dynamics for EUROFER-97 under irradiation" \
      org.opencontainers.image.source="https://github.com/Ghoniem/RadCluster" \
      org.opencontainers.image.licenses="MIT" \
      org.opencontainers.image.authors="Nasr M. Ghoniem"

ENV DEBIAN_FRONTEND=noninteractive
# libsundials-dev rather than the versioned runtime packages: their names carry
# the Ubuntu time_t transition suffix (libsundials-cvode6t64 and friends) and
# shift between releases, which would break this build on the next base-image
# bump for the sake of a few MB.  The -dev package depends on the runtime ones,
# so this is correct by construction.
RUN apt-get update && apt-get install -y --no-install-recommends \
        python3 \
        python3-pip \
        libsundials-dev \
        liblapack3 \
        libopenblas0 \
        libsuitesparse-dev \
    && rm -rf /var/lib/apt/lists/*

# Ubuntu 24.04 marks the system interpreter externally managed (PEP 668).  This
# image is a single-purpose appliance, not a shared system, so installing into
# it is the intent rather than an accident.
COPY requirements.txt /tmp/requirements.txt
RUN pip3 install --no-cache-dir --break-system-packages -r /tmp/requirements.txt \
    && rm /tmp/requirements.txt

WORKDIR /work
COPY --from=build /out/solver /usr/local/bin/solver
COPY radcluster_code/py_utils   ./py_utils
COPY radcluster_code/input      ./input
COPY radcluster_code/codes      ./codes
COPY README.md LICENSE          ./

# cpp_bridge.py looks for the compiled binary under build/; point it at the one
# installed on PATH rather than copying it twice.
RUN mkdir -p build/Release && ln -s /usr/local/bin/solver build/Release/solver

ENV PYTHONPATH=/work \
    PYTHONIOENCODING=utf-8 \
    MPLBACKEND=Agg

# Fails the build if the RAG cannot be declared from the shipped workbook, so a
# broken image is never published.
RUN python3 codes/Python_Testing/check_eurofer_rag.py > /dev/null \
    && echo "RAG declaration OK"

CMD ["/bin/bash"]
