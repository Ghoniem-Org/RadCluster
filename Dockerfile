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
        gfortran \
        cmake \
        ca-certificates \
        curl \
        liblapack-dev \
        libopenblas-dev \
        libsuitesparse-dev \
    && rm -rf /var/lib/apt/lists/*

# SUNDIALS from source, not from apt.  Ubuntu 24.04 ships 6.4.1, and this
# solver is written against 7: core/solver.cpp calls
#     SUNContext_Create(SUN_COMM_NULL, &sunctx)
# and SUN_COMM_NULL does not exist before v7 (v6 took a raw void* comm).
# Building against the distribution package fails to compile, and pinning
# 7.1.1 here matches the version CLAUDE.md pins for every other platform
# rather than quietly running the image on a different SUNDIALS.
# gfortran is in that list for SUNDIALS, not for this project: enabling LAPACK
# makes SUNDIALS probe the Fortran name-mangling scheme by compiling a small
# Fortran library (cmake/tpl/SundialsLapack.cmake), and without a Fortran
# compiler that probe fails and the configure aborts with
#     FATAL_ERROR: SUNDIALS interface to LAPACK is not functional.
# which is what broke the first two attempts at this image.

ARG SUNDIALS_VERSION=7.1.1
RUN curl -fsSL -o /tmp/sundials.tar.gz \
        "https://github.com/LLNL/sundials/releases/download/v${SUNDIALS_VERSION}/sundials-${SUNDIALS_VERSION}.tar.gz" \
    && tar -xzf /tmp/sundials.tar.gz -C /tmp \
    && cmake -S "/tmp/sundials-${SUNDIALS_VERSION}" -B /tmp/sundials-build \
        -DCMAKE_BUILD_TYPE=Release \
        -DCMAKE_INSTALL_PREFIX=/opt/sundials \
        -DBUILD_SHARED_LIBS=ON \
        -DBUILD_STATIC_LIBS=OFF \
        -DEXAMPLES_ENABLE_C=OFF \
        -DEXAMPLES_INSTALL=OFF \
        -DENABLE_LAPACK=ON \
        -DENABLE_KLU=ON \
        -DKLU_INCLUDE_DIR=/usr/include/suitesparse \
        -DKLU_LIBRARY_DIR=/usr/lib/x86_64-linux-gnu \
    && cmake --build /tmp/sundials-build --parallel "$(nproc)" \
    && cmake --install /tmp/sundials-build \
    && rm -rf /tmp/sundials.tar.gz "/tmp/sundials-${SUNDIALS_VERSION}" /tmp/sundials-build

WORKDIR /src
COPY radcluster_code/cpp_utils ./cpp_utils

# -march=native would bake this builder's instruction set into a binary that
# other people pull and run; x86-64-v2 (SSE4.2 / POPCNT, 2009 and later) is the
# portable floor.  Fast-math stays off here for the same reason it is off in the
# local build: it perturbs IEEE-754 semantics and can make CVODE's error
# estimator diverge on these stiff systems.
RUN cmake -S cpp_utils -B build \
        -DCMAKE_BUILD_TYPE=Release \
        -DCMAKE_PREFIX_PATH=/opt/sundials \
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
RUN apt-get update && apt-get install -y --no-install-recommends \
        python3 \
        python3-pip \
        liblapack3 \
        libopenblas0 \
        libsuitesparse-dev \
    && rm -rf /var/lib/apt/lists/*

# The SUNDIALS shared libraries come from the build stage, not from apt, for
# the version reason above -- apt would install 6.4.1 beside a binary linked
# against 7.1.1.
COPY --from=build /opt/sundials/lib /opt/sundials/lib
ENV LD_LIBRARY_PATH=/opt/sundials/lib

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
# The link gate is the one that would have caught the SUNDIALS version
# mismatch at image-build time rather than at someone's first run.
RUN if ldd /usr/local/bin/solver | grep -q "not found"; then \
        echo "solver has unresolved shared libraries:"; \
        ldd /usr/local/bin/solver; exit 1; \
    else echo "solver links OK"; fi
RUN python3 codes/Python_Testing/check_eurofer_rag.py > /dev/null \
    && echo "RAG declaration OK"

CMD ["/bin/bash"]
