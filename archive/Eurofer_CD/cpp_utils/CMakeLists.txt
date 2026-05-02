cmake_minimum_required(VERSION 3.10)
project(eurofer_cd_solver)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# ── SUNDIALS ──────────────────────────────────────────────────────────────────
# Same Libraries/ layout as Full_CD:
#   <repo_root>/Libraries/sundials-7.1.1/
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

# ── Executable ────────────────────────────────────────────────────────────────
# solver.cpp lives in Eurofer_CD/codes/
# rate_equations.cpp lives here in cpp_utils/
add_executable(solver
    ../codes/solver.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/rate_equations.cpp
)

# Include cpp_utils/ so solver.cpp finds parameters.h and rate_equations.h
target_include_directories(solver PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

# ── Compiler optimizations ────────────────────────────────────────────────────
# NOTE: -ffast-math / /fp:fast are intentionally omitted — they alter IEEE-754
# floating-point semantics and can cause CVODE's error estimator to diverge on
# stiff problems.
if(MSVC)
    target_compile_options(solver PRIVATE /arch:AVX2 /GL /fp:precise)
    target_link_options(solver PRIVATE /LTCG)
else()
    target_compile_options(solver PRIVATE -O3 -march=native)
endif()

# ── OpenMP (Phase IV: Multithread-OpenMP) ─────────────────────────────────────
if(APPLE AND CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    execute_process(
        COMMAND brew --prefix libomp
        OUTPUT_VARIABLE _BREW_LIBOMP
        OUTPUT_STRIP_TRAILING_WHITESPACE
        ERROR_QUIET)
    if(_BREW_LIBOMP AND EXISTS "${_BREW_LIBOMP}")
        set(OpenMP_CXX_FLAGS          "-Xpreprocessor -fopenmp")
        set(OpenMP_CXX_LIB_NAMES      "omp")
        set(OpenMP_omp_LIBRARY        "${_BREW_LIBOMP}/lib/libomp.dylib")
        include_directories(SYSTEM "${_BREW_LIBOMP}/include")
    endif()
endif()

find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(solver OpenMP::OpenMP_CXX)
    target_compile_definitions(solver PRIVATE CD_HAVE_OPENMP=1)
    message(STATUS "OpenMP found — Phase IV (Multithread-OpenMP) ENABLED")
else()
    message(WARNING
        "OpenMP not found — Phase IV disabled (window_mode=4 will fall back to Phase III).\n"
        "  On macOS:   brew install libomp  then re-run cmake.\n"
        "  On Windows/MSVC: OpenMP is bundled — this should not happen.")
endif()

# ── Link SUNDIALS — CVODE + ARKODE + linear solvers ──────────────────────────
target_link_libraries(solver
    SUNDIALS::cvode
    SUNDIALS::arkode
    SUNDIALS::nvecserial
    SUNDIALS::sunlinsolband
    SUNDIALS::sunlinsolspgmr
)
