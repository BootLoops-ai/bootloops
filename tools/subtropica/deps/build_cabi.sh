#!/bin/bash
# build_cabi.sh — build of the patched, static-isolated HyperFLINT C-ABI lib
# for the subtropica in-process (ccall) transport (see deps/BUILD.md).
# Source tree: the patched work tree (pristine clone + the F2 patch ONLY).
set -e
ulimit -v 32505856
DEPS="$(cd "$(dirname "$0")" && pwd)"
# HF_CABI_WORK must point at a work tree holding: fix/HyperFLINT-f2 (pristine
# HyperFLINT clone + the F2 patch ONLY), hf_cabi.map, flint-static-install/,
# deps-install/ (static PIC FLINT/GMP/MPFR). Refuse loudly if unset.
: "${HF_CABI_WORK:?set HF_CABI_WORK to your HyperFLINT C-ABI work tree (see GUIDE.md)}"
CABI=$HF_CABI_WORK
SRC=$CABI/fix/HyperFLINT-f2
# pkg-config path for a static PIC FLINT install
: "${HF_FLINT_PKGCONFIG:?set HF_FLINT_PKGCONFIG to your flint-install/lib/pkgconfig}"
export PKG_CONFIG_PATH=$HF_FLINT_PKGCONFIG
# -DHF_VERSION=1.2.8 MANDATORY: source builds stamp 1.2.0 by default;
# the dlopen gate checks strict equality with 1.2.8.
cmake -S "$SRC" -B "$DEPS/build" -DCMAKE_BUILD_TYPE=Release \
      -DHF_VERSION=1.2.8 -DHF_LIBRARYLINK=OFF -DHF_MIMALLOC=OFF \
      -DCMAKE_POSITION_INDEPENDENT_CODE=ON > "$DEPS/build.cfg.log" 2>&1
cmake --build "$DEPS/build" --target hyperflint -j 24 > "$DEPS/build.log" 2>&1
# static-isolated link: PRIVATE static PIC FLINT/GMP/MPFR (never Julia's
# GC-hooked libgmp — see deps/BUILD.md); exactly 8 exported symbols.
SYMS="-u hf_free_string -u hf_partial_fractions -u hf_linear_factors \
      -u hf_find_lr_orders -u hf_find_lr_orders_scan -u hf_factor_table \
      -u hf_hyperflint_sym -u hf_version_string"
cp "$CABI/hf_cabi.map" "$DEPS/hf_cabi.map"
g++ -shared -o "$DEPS/libhf_cabi_static_f2.so" \
    -Wl,--version-script="$DEPS/hf_cabi.map" $SYMS \
    "$DEPS/build/libhyperflint.a" \
    "$CABI/flint-static-install/lib/libflint.a" \
    "$CABI/deps-install/lib/libmpfr.a" \
    "$CABI/deps-install/lib/libgmp.a" \
    -fopenmp -lpthread -lm
echo BUILD_OK
