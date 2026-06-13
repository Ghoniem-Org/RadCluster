cmake_minimum_required(VERSION 3.10)
project(radcluster_2_1_solver)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# ══════════════════════════════════════════════════════════════════════════════
# RadCluster_2_1 C++ solver — two-layer core/materials architecture.
#
#   cpp_utils/core/                 host-independent solver layer
#     solver.cpp          — CVODE driver / main()
#     rhs_dispatch.cpp    — preconditioner dispatch + Jacobi + Woodbury
#     sparse_jacobian.cpp — colored-FD sparse Jacobian (KLU path)
#     parameters.h        — parameter parsing
#     rate_equations.h    — RHS-dispatch scaffolding (shared header)
#     sparse_jacobian.h
#
#   cpp_utils/materials/eurofer97/  EUROFER97-specific physics layer
#     rate_kernels.cpp    — P1-P8 rate-equation kernel arithmetic
#                           (rhs_full_CD / rhs_case1 / rhs_case2 / rhs_bin_moment)
#
# This mirrors the Python py_utils core/materials split.  Swapping the material
# model means adding a new materials/<name>/ directory; the core layer is
# untouched.
# ══════════════════════════════════════════════════════════════════════════════

set(CORE_DIR      "${CMAKE_CURRENT_SOURCE_DIR}/core")
set(MATERIALS_DIR "${CMAKE_CURRENT_SOURCE_DIR}/materials")
set(EUROFER_DIR   "${MATERIALS_DIR}/eurofer97")

# ── SUNDIALS 7.1.1 ────────────────────────────────────────────────────────────
# Mirrors the Eurofer_CD CMakeLists layout: Libraries/ sibling of repository root.
get_filename_component(LIBRARIES_DIR
    "${CMAKE_CURRENT_SOURCE_DIR}/../../../Libraries" ABSOLUTE)
if(APPLE)
    list(PREPEND CMAKE_PREFIX_PATH "${LIBRARIES_DIR}/sundials-7.1.1/mac_install")
elseif(WIN32)
    list(PREPEND CMAKE_PREFIX_PATH
        "${LIBRARIES_DIR}/sundials-7.1.1/win_install"
        "${LIBRARIES_DIR}/sundials-7.1.1")
else()
    list(PREPEND CMAKE_PREFIX_PATH "${LIBRARIES_DIR}/sundials-7.1.1")
endif()
find_package(SUNDIALS REQUIRED)

# ── Solver executable ─────────────────────────────────────────────────────────
add_executable(solver
    ${CORE_DIR}/solver.cpp
    ${CORE_DIR}/rhs_dispatch.cpp
    ${CORE_DIR}/sparse_jacobian.cpp
    ${EUROFER_DIR}/rate_kernels.cpp
)

# core/ carries every shared header (parameters.h, rate_equations.h,
# sparse_jacobian.h).  materials/ is on the path too so future material
# models can ship their own headers.
target_include_directories(solver PRIVATE ${CORE_DIR} ${MATERIALS_DIR})

# ── Compiler optimizations ────────────────────────────────────────────────────
# NOTE: -ffast-math / /fp:fast are intentionally omitted — they alter IEEE-754
# semantics and can cause CVODE's error estimator to diverge on stiff problems.
if(MSVC)
    target_compile_options(solver PRIVATE /arch:AVX2 /GL /fp:precise)
    target_link_options(solver PRIVATE /LTCG)
else()
    target_compile_options(solver PRIVATE -O3 -march=native)
endif()

# ── OpenMP (active_window solver mode) ────────────────────────────────────────
if(APPLE AND CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    execute_process(
        COMMAND brew --prefix libomp
        OUTPUT_VARIABLE _BREW_LIBOMP
        OUTPUT_STRIP_TRAILING_WHITESPACE
        ERROR_QUIET)
    if(_BREW_LIBOMP AND EXISTS "${_BREW_LIBOMP}")
        set(OpenMP_CXX_FLAGS     "-Xpreprocessor -fopenmp")
        set(OpenMP_CXX_LIB_NAMES "omp")
        set(OpenMP_omp_LIBRARY   "${_BREW_LIBOMP}/lib/libomp.dylib")
        include_directories(SYSTEM "${_BREW_LIBOMP}/include")
    endif()
endif()

find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(solver OpenMP::OpenMP_CXX)
    target_compile_definitions(solver PRIVATE CD_HAVE_OPENMP=1)
    message(STATUS "OpenMP found — active_window mode ENABLED")
else()
    message(WARNING
        "OpenMP not found — active_window (window_mode=4) will run serial"
        " (single-threaded sliding window).\n"
        "  macOS:   brew install libomp  then re-run cmake.\n"
        "  Windows: OpenMP is bundled with MSVC.")
endif()

# ── LAPACK (Woodbury preconditioner: dgbtrf/dgbtrs, dgetrf/dgetrs) ───────────
# RadCluster_1_0 located LAPACK via the vcpkg toolchain
# (C:/vcpkg/installed/x64-windows/lib/{lapack,openblas}.lib).  RadCluster_2_1
# previously built WITHOUT the toolchain flag, so find_package(LAPACK) failed
# and the Woodbury preconditioner was compiled out.
#
# This block reproduces the v1.0 result without depending on the caller
# remembering -DCMAKE_TOOLCHAIN_FILE: it seeds the vcpkg install tree onto the
# search paths, runs find_package(LAPACK), and — if that still fails — links
# the vcpkg import libraries directly.
set(LAPACK_FOUND FALSE)
if(APPLE)
    # macOS Accelerate framework includes LAPACK
    find_library(ACCELERATE_FRAMEWORK Accelerate)
    if(ACCELERATE_FRAMEWORK)
        message(STATUS "Using macOS Accelerate framework for LAPACK")
    endif()
else()
    # Seed common vcpkg locations so find_package(LAPACK) succeeds even when
    # the vcpkg toolchain file was not passed at configure time.
    set(_VCPKG_CANDIDATES
        "$ENV{VCPKG_ROOT}/installed/x64-windows"
        "C:/vcpkg/installed/x64-windows")
    foreach(_vc ${_VCPKG_CANDIDATES})
        if(EXISTS "${_vc}/lib")
            list(APPEND CMAKE_PREFIX_PATH  "${_vc}")
            list(APPEND CMAKE_LIBRARY_PATH "${_vc}/lib")
            message(STATUS "LAPACK search: added vcpkg prefix ${_vc}")
        endif()
    endforeach()
    find_package(LAPACK)
endif()

if(APPLE AND ACCELERATE_FRAMEWORK)
    target_link_libraries(solver ${ACCELERATE_FRAMEWORK})
    target_compile_definitions(solver PRIVATE CD_HAVE_LAPACK=1)
    message(STATUS "LAPACK (Accelerate) — Woodbury preconditioner ENABLED")
elseif(LAPACK_FOUND)
    target_link_libraries(solver ${LAPACK_LIBRARIES})
    target_compile_definitions(solver PRIVATE CD_HAVE_LAPACK=1)
    message(STATUS "LAPACK found — Woodbury preconditioner ENABLED")
else()
    # find_package(LAPACK) did not succeed — fall back to the vcpkg import
    # libraries directly (the exact pair RadCluster_1_0 linked against).
    set(_LAPACK_LIB "")
    set(_BLAS_LIB   "")
    foreach(_vc ${_VCPKG_CANDIDATES})
        if(NOT _LAPACK_LIB AND EXISTS "${_vc}/lib/lapack.lib")
            set(_LAPACK_LIB "${_vc}/lib/lapack.lib")
        endif()
        if(NOT _BLAS_LIB AND EXISTS "${_vc}/lib/openblas.lib")
            set(_BLAS_LIB "${_vc}/lib/openblas.lib")
        endif()
    endforeach()
    if(_LAPACK_LIB AND _BLAS_LIB)
        target_link_libraries(solver "${_LAPACK_LIB}" "${_BLAS_LIB}")
        target_compile_definitions(solver PRIVATE CD_HAVE_LAPACK=1)
        message(STATUS
            "LAPACK found via vcpkg import libs — Woodbury preconditioner ENABLED\n"
            "  LAPACK: ${_LAPACK_LIB}\n"
            "  BLAS:   ${_BLAS_LIB}")
    else()
        message(WARNING
            "LAPACK not found — Woodbury preconditioner will be COMPILED OUT.\n"
            "  Solver will fall back to Jacobi preconditioner at runtime if\n"
            "  prec_type=1 is requested.\n"
            "  Windows: install LAPACK via vcpkg ('vcpkg install lapack-reference'\n"
            "    or 'vcpkg install openblas[lapack]'); the build auto-detects\n"
            "    C:/vcpkg/installed/x64-windows or $ENV{VCPKG_ROOT}.")
    endif()
endif()

# ── Link SUNDIALS ─────────────────────────────────────────────────────────────
target_link_libraries(solver
    SUNDIALS::cvode
    SUNDIALS::nvecserial
    SUNDIALS::sunlinsolband
    SUNDIALS::sunlinsolspgmr
    SUNDIALS::sunlinsoldense
)

# ── KLU sparse direct solver (linsol=3) ──────────────────────────────────────
# SUNDIALS::sunlinsolklu is only present if SUNDIALS was built with
# -DENABLE_KLU=ON and SuiteSparse/KLU is available. If absent, the solver
# still builds and the linsol=3 branch errors out at runtime with a
# helpful message.
if(TARGET SUNDIALS::sunlinsolklu)
    target_link_libraries(solver SUNDIALS::sunlinsolklu)
    # The KLU path uses SUNSparseMatrix, which lives in sunmatrixsparse.
    if(TARGET SUNDIALS::sunmatrixsparse)
        target_link_libraries(solver SUNDIALS::sunmatrixsparse)
    endif()
    target_compile_definitions(solver PRIVATE CD_HAVE_KLU=1)
    message(STATUS "KLU found — linsol='klu' (sparse direct) ENABLED")
else()
    message(WARNING
        "KLU not found — linsol='klu' (sparse direct solver) DISABLED.\n"
        "  To enable: rebuild SUNDIALS with -DENABLE_KLU=ON and ensure\n"
        "  SuiteSparse/KLU is available (Linux/Mac: libsuitesparse-dev or\n"
        "  brew install suite-sparse; Windows: vcpkg install suitesparse).\n"
        "  Then re-run cmake on this project.")
endif()

# ── Windows runtime-DLL deployment (POST_BUILD) ──────────────────────────────
# The solver links SUNDIALS (DYNAMIC, from Libraries/sundials-7.1.1/win_install)
# and the vcpkg LAPACK/OpenBLAS + KLU/SuiteSparse import libs.  None of those
# DLLs sit on PATH, so a fresh build — or the notebook's build cell, which does
# NOT pass -DCMAKE_TOOLCHAIN_FILE and so skips vcpkg's applocal deployment —
# produces a solver.exe that fails to load with 0xC0000135
# (STATUS_DLL_NOT_FOUND).  Copy every required runtime DLL next to the exe after
# each build so the binary is self-contained on a clean checkout.
#
# Resolved at CONFIGURE time via file(GLOB); missing source trees are skipped,
# so this is a no-op (not an error) on a machine where a path is absent.
if(WIN32)
    set(_RUNTIME_DLLS "")

    # (1) SUNDIALS shared libraries — whole bin/ of the install tree we linked.
    set(_SUNDIALS_BIN "${LIBRARIES_DIR}/sundials-7.1.1/win_install/bin")
    if(EXISTS "${_SUNDIALS_BIN}")
        file(GLOB _SUN_DLLS "${_SUNDIALS_BIN}/*.dll")
        list(APPEND _RUNTIME_DLLS ${_SUN_DLLS})
    endif()

    # (2) vcpkg runtime DLLs — OpenBLAS/LAPACK, SuiteSparse/KLU, and the MinGW
    # runtime that OpenBLAS depends on.  Curated patterns (NOT the whole vcpkg
    # bin).  First existing vcpkg tree from the LAPACK search above wins.
    foreach(_vc ${_VCPKG_CANDIDATES})
        if(EXISTS "${_vc}/bin")
            file(GLOB _VCPKG_DLLS
                "${_vc}/bin/openblas.dll"
                "${_vc}/bin/lapack.dll"
                "${_vc}/bin/liblapack.dll"
                "${_vc}/bin/klu.dll"
                "${_vc}/bin/amd.dll"
                "${_vc}/bin/btf.dll"
                "${_vc}/bin/colamd.dll"
                "${_vc}/bin/suitesparseconfig.dll"
                "${_vc}/bin/libgcc_s_seh-1.dll"
                "${_vc}/bin/libgfortran-5.dll"
                "${_vc}/bin/libquadmath-0.dll"
                "${_vc}/bin/libwinpthread-1.dll")
            list(APPEND _RUNTIME_DLLS ${_VCPKG_DLLS})
            break()
        endif()
    endforeach()

    list(REMOVE_DUPLICATES _RUNTIME_DLLS)
    if(_RUNTIME_DLLS)
        list(LENGTH _RUNTIME_DLLS _N_DLLS)
        add_custom_command(TARGET solver POST_BUILD
            COMMAND ${CMAKE_COMMAND} -E copy_if_different
                    ${_RUNTIME_DLLS} "$<TARGET_FILE_DIR:solver>"
            COMMAND_EXPAND_LISTS
            COMMENT "Deploying runtime DLLs next to solver.exe")
        message(STATUS
            "Runtime DLL deploy: ${_N_DLLS} DLLs will be copied next to solver.exe")
    else()
        message(WARNING
            "Runtime DLL deploy: no SUNDIALS/vcpkg DLLs found to copy.\n"
            "  solver.exe may fail to load at runtime (0xC0000135) unless the\n"
            "  DLL directories are on PATH.")
    endif()
endif()
