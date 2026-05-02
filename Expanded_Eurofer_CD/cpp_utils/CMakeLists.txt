cmake_minimum_required(VERSION 3.10)
project(expanded_eurofer_cd_solver)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

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
# solver.cpp lives in codes/
# rate_equations.cpp lives here in cpp_utils/
add_executable(solver
    ../codes/solver.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/rate_equations.cpp
)

target_include_directories(solver PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

# ── Compiler optimizations ────────────────────────────────────────────────────
# NOTE: -ffast-math / /fp:fast are intentionally omitted — they alter IEEE-754
# semantics and can cause CVODE's error estimator to diverge on stiff problems.
if(MSVC)
    target_compile_options(solver PRIVATE /arch:AVX2 /GL /fp:precise)
    target_link_options(solver PRIVATE /LTCG)
else()
    target_compile_options(solver PRIVATE -O3 -march=native)
endif()

# ── OpenMP (sliding_OpenMP solver mode) ───────────────────────────────────────
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
    message(STATUS "OpenMP found — sliding_OpenMP mode ENABLED")
else()
    message(WARNING
        "OpenMP not found — sliding_OpenMP (window_mode=4) will fall back to"
        " cpp_sliding_win.\n"
        "  macOS:   brew install libomp  then re-run cmake.\n"
        "  Windows: OpenMP is bundled with MSVC.")
endif()

# ── LAPACK (Woodbury preconditioner: dgbtrf/dgbtrs, dgetrf/dgetrs) ───────────
if(APPLE)
    # macOS Accelerate framework includes LAPACK
    find_library(ACCELERATE_FRAMEWORK Accelerate)
    if(ACCELERATE_FRAMEWORK)
        message(STATUS "Using macOS Accelerate framework for LAPACK")
    endif()
else()
    find_package(LAPACK)
endif()

# ── Link SUNDIALS + LAPACK ───────────────────────────────────────────────────
target_link_libraries(solver
    SUNDIALS::cvode
    SUNDIALS::arkode
    SUNDIALS::nvecserial
    SUNDIALS::sunlinsolband
    SUNDIALS::sunlinsolspgmr
    SUNDIALS::sunlinsoldense
)
if(APPLE AND ACCELERATE_FRAMEWORK)
    target_link_libraries(solver ${ACCELERATE_FRAMEWORK})
    target_compile_definitions(solver PRIVATE CD_HAVE_LAPACK=1)
elseif(LAPACK_FOUND)
    target_link_libraries(solver ${LAPACK_LIBRARIES})
    target_compile_definitions(solver PRIVATE CD_HAVE_LAPACK=1)
else()
    message(WARNING
        "LAPACK not found — Woodbury preconditioner will be COMPILED OUT.\n"
        "  Solver will fall back to Jacobi preconditioner at runtime if\n"
        "  prec_type=1 is requested.\n"
        "  Windows: install LAPACK via vcpkg ('vcpkg install lapack-reference'\n"
        "    or 'vcpkg install openblas[lapack]') and pass\n"
        "    -DCMAKE_TOOLCHAIN_FILE=<vcpkg>/scripts/buildsystems/vcpkg.cmake\n"
        "    to cmake, or install Intel MKL.")
endif()
