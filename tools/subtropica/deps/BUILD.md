# deps/ — the patched HyperFLINT C-ABI library (in-process transport)

## Artifact

`libhf_cabi_static_f2.so` — a statically isolated shared library exposing
HyperFLINT's 8-symbol stable C ABI, built from a **patched** upstream
source tree: the pristine HyperFLINT clone at the pinned commit plus one
patch, referred to throughout this package as the F2 patch — a per-call
PolyCtx hygiene fix: we observed heap corruption when reusing the v1.2.8
library in-process in our build; the CLI transport (one process per call)
is unaffected and is the supported route. The patch is not distributed with this release; without a
patched tree the in-process transport is unavailable and the CLI transport
(the default) is the supported route. The 56 MB prebuilt binary is not
included either.

- Exactly 8 exported symbols (version script `hf_cabi.map`): hf_free_string,
  hf_partial_fractions, hf_linear_factors, hf_find_lr_orders,
  hf_find_lr_orders_scan, hf_factor_table, hf_hyperflint_sym,
  hf_version_string. Everything else hidden — Nemo/FLINT_jll safe. Ops
  without an ABI entry point (parse_expr, eval, the period and Vieta ops)
  are the "C-ABI gap G1" named in `src/hf_bridge.jl`: they always go
  through the CLI, even under `transport=:cabi`.
- Ownership: every response `char*` is freed by the caller via
  `hf_free_string`; a NULL return is a catastrophic allocation failure.
- PRIVATE static PIC FLINT 3.5.0 + GMP 6.3.0 + MPFR 4.2.1 linked in
  (archives from the work tree's `flint-static-install`/`deps-install`).
  NEVER Julia's GC-hooked libgmp: Julia installs GMP memory hooks, so a
  library sharing Julia's libgmp corrupts memory across the boundary.
  Dynamic deps: libstdc++/libgomp/libm/libc only (check with ldd).
- `-DHF_VERSION=1.2.8` at cmake is MANDATORY (source builds stamp 1.2.0 by
  default; the dlopen gate in `src/hf_bridge.jl` checks strict equality).

## Rebuild

`./build_cabi.sh` — out-of-tree cmake, then a static-isolated link; needs
`HF_CABI_WORK` (a work tree holding `fix/HyperFLINT-f2` = the patched
clone, `hf_cabi.map`, `flint-static-install/`, `deps-install/`) and
`HF_FLINT_PKGCONFIG`; the script refuses loudly when either is unset. The
cmake `build/` tree may be deleted after linking — only the .so and the
logs need to stay here.

## Post-build acceptance (rerun after ANY rebuild)

`ulimit -v 32505856; julia deps/verify_cabi.jl` — all checks must PASS:
version string == 1.2.8 strict, envelope `"hf_version":"1.2.8"`,
unary-minus canary through the in-process parser, the F2 pair-repro
sequence error-free with byte-identical repeats (pfrac AND linear_factors
faces). An UNPATCHED library fails or corrupts on the pair sequence
(heap-layout dependent — an ASAN build makes it deterministic).

## Use

Consumed by `src/hf_bridge.jl` `transport=:cabi` (OPT-IN; the CLI stays
the default and the validation oracle). Override the path via
`SUBTROPICA_HF_CABI_LIB`. RULES: one dlopen handle per PROCESS; ccalls
serialized (threads on a shared handle crash inside FLINT — measured);
farm parallelism = worker processes, ended by PID; timeouts are
per-process, not per-call.
