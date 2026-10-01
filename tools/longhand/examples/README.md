# Longhand — worked example of route B: `ndpent_top`

`ndpent_top/` holds the records of one complete run of the disteval chain
(`../disteval/pysecdec_point.py`, manual `../DISTEVAL.md`) on the generic-mass five-point
top-sector scalar `ndpent[1,1,1,1,1,1,1,1,0,0,0]` — the double pentagon with all eight
propagators, 332 sectors, 968 kernels after pySecDec 1.6.6's geometric decomposition — at
the Euclidean point `s12=-7 s23=-9 s34=-13 s45=-10 s15=-5 mm3=-12 mm4=-8`:

    eps^0 = 8.246341248953e-01 +- 7.3e-08     (stage B; per-order digits 9/7/7/7/7 for eps^-4..eps^0; digit class 7)

Every file is sha256-pinned in `ndpent_top/PINS.json`, and the directory holds nothing
else (the battery asserts both). What each record is:

| file | what it records |
|---|---|
| `RESULT_A.json`, `RESULT_B.json`, `RESULT_C.json` | the three assembled stage results (lattice 1e4; lattice 1e5 in 28 chunks; the median-lattice second integrator in 28 chunks), per regulator order with bars |
| `COMPARE_record.json` | the pairwise per-order compare of the three |
| `GATE_record.json` | the acceptance check over them: the file `--stage check --example` must reproduce (21 checks by name) |
| `CHUNKS_B_record.json`, `SETTINGS_B_record.json`, `PROGRESS_record.jsonl` | stage B's chunk manifest, settings (per-chunk epsabs `1.706e-10`) and per-chunk progress lines (walls, maxrss, disteval statistics) |
| `smoke_build_host.json`, `smoke_second_host.json`, `SMOKE_build_host_receipt.json` | the two coarse independent evaluations of the same library that serve as the FOREIGN control, and the build host's receipt |
| `BASIS_record.json`, `FENCE_AT_SPAWN_record.json`, `LEAF_STAMP_record.txt` | the measured memory basis (per-worker and driver maxrss), the cgroup readback at spawn, and the end-of-run stamp the check reads (`oom_kill 0`, `memory.max 377241600`) |
| `GENERATE_RECEIPT_record.json` | the `loop_package` generate step's receipt |
| `package_ndpent_top.json`, `package_ndpent_top_integral.json` | the compiled package's two descriptor JSONs (parameter names, prefactor expansion, sector/order kernel lists — pySecDec-generated metadata for our integral, no pySecDec code), so the chunk split, the stand-in stage and the pins legs run without the compiled package |

The compiled package itself (the `.so` libraries and `coefficients/`, 13 MB, a build of
its host) is not shipped; `../disteval/PACKAGE_PINS.json` pins its five files so a rebuilt
copy can be verified (`pysecdec_point.py --check-package --disteval-dir PKG/disteval`).
To run your own family the same way, see "Your own family" in `../DISTEVAL.md`.
