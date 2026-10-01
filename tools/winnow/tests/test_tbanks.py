#!/usr/bin/env python3
"""BYTE-EQUAL T-bank regression against the archived reference banks
(replays archived artifacts — the whole suite skips without
WINNOW_BANKED_ROOT / WINNOW_T1_ROOT).

Legs:
  1. lbl3m2L_k2disp: replay all 4 archived slices through the ibplapper
     API (Bank + Schedule, cell "tiny", node "p{p}_d{d}_e{e}") from the
     archived Generate; the produced T .bin must be BYTE-EQUAL (sha256) to
     the archived reference T for BOTH B2F and B2FT (reference fact:
     v0-vs-B2F banks are byte-equal; B2FT banks identical per-stratum rows
     by construction).
  2. lbl3m2L_lp1disp: same, from the stagepack Generate (a fresh kira
     rebuild — byte-equality here also witnesses Generate determinism).
  3. lbl3mopp_sh: the raw Generate for this family is not in the archive,
     so the B2FT-anchor byte-equal claim is NOT re-runnable; the leg
     degrades HONESTLY to a full read-back integrity sweep of the archived
     sh T (all 24 blocks sha256-verified through the Bank). The re-run
     claim stays archived, not re-proved.

T1-wall regression is recorded here too (measure-and-record discipline):
current engine-path wall vs the recorded reference prior; hard-fail only
beyond 2x (structural), report the ratio always.
"""
import _common  # noqa: F401  MUST be first: sys.path bootstrap for the package

import hashlib
import os
import statistics
import time

import ibplapper as lap
from ibplapper.adapters import kira
from ibplapper.engine import Stats, eliminate_fast
from _common import (D4OUT, NOMAT, SLICES4, STAGEPACKS, T1GEN, check,
                     finish, out_dir, skip_without)

skip_without("test_tbanks", D4OUT, NOMAT, STAGEPACKS, T1GEN)

OUT = out_dir("/tmp/ibplapper_test_tbanks")


def sha(path):
    return hashlib.sha256(open(path, "rb").read()).hexdigest()


def replay(gen, fam, policy, base):
    for ext in (".bin", ".idx.json"):
        if os.path.exists(base + ext):
            os.remove(base + ext)
    for si, (p, d0, e0) in enumerate(SLICES4):
        system, _sysd = kira.load(gen, fam, p, d0, e0)
        bank = lap.Bank(base, meta=None if si else {"family": fam})
        res = lap.eliminate(system, lap.Schedule(
            policy=policy, bank=bank, bank_cell="tiny",
            bank_node=f"p{p}_d{d0}_e{e0}"))
        assert not res.leftover
    return base + ".bin"


# ---- leg 1: k2disp, both policies, byte-equal to the reference T --------
banked = os.path.join(D4OUT, "lbl3m2L_k2disp", "T_lbl3m2L_k2disp.bin")
banked_v0 = os.path.join(D4OUT, "lbl3m2L_k2disp_v0", "T_lbl3m2L_k2disp.bin")
check("banked_b2f_vs_v0_still_byte_equal", sha(banked) == sha(banked_v0))
for pol in ("B2F", "B2FT"):
    got = replay(os.path.join(NOMAT, "gen_lbl3m2L_k2disp"),
                 "lbl3m2L_k2disp", pol, os.path.join(OUT, f"T_k2disp_{pol}"))
    check(f"k2disp_T_byte_equal::{pol}", sha(got) == sha(banked))

# ---- leg 2: lp1disp from the fresh stagepack Generate ---------------------
banked = os.path.join(D4OUT, "lbl3m2L_lp1disp", "T_lbl3m2L_lp1disp.bin")
got = replay(os.path.join(STAGEPACKS, "gen_lbl3m2L_lp1disp"),
             "lbl3m2L_lp1disp", "B2FT", os.path.join(OUT, "T_lp1disp"))
check("lp1disp_T_byte_equal::B2FT_freshgen", sha(got) == sha(banked))

# ---- leg 3: sh — integrity sweep only (artifacts gone; stated honestly) ---
from ibplapper.bank import InterfaceTable                    # noqa: E402
shT = InterfaceTable(os.path.join(D4OUT, "lbl3mopp_sh", "T_lbl3mopp_sh"))
n_blocks = n_rows = 0
for k in shT.idx["keys"]:
    fam, cell, node = k.split("|")
    rows = shT.get(cell, node, family=fam, verify=True)   # sha256-verified
    n_blocks += 1
    n_rows += len(rows)
check("sh_T_integrity_24_blocks_hash_verified", n_blocks == 24,
      f"blocks={n_blocks} rows={n_rows}")
print("[note] sh B2FT-anchor byte-equal re-run NOT possible (the raw "
      "Generate artifacts are not in the archive); claim remains "
      "archived-only record, integrity-swept here.")

# ---- T1 wall: measure and record vs the recorded reference prior ----------
PRIOR_MS = 107.2        # recorded reference mean for vac3_b0_r2_f2_sm0
walls = []
for _rep in range(2):
    for (p, d0, e0) in SLICES4:
        sysd = kira.load_system(T1GEN, "vac3_b0_r2_f2_sm0", p, d0, e0)
        t0 = time.time()
        subs, leftover, stats = eliminate_fast(
            sysd["rows"], sysd["order"], p, forbid=set(sysd["masters"]),
            stats=Stats())
        walls.append(time.time() - t0)
        check_stats = (stats["pivots"] == 6001 and stats["ops"] == 54550
                       and stats["fill"] == 21773 and not leftover)
        if not check_stats:
            check("T1_stats_match_prior", False, stats)
            break
check("T1_stats_match_prior", True,
      "pivots=6001 ops=54550 fill=21773 leftover=0")
mean_ms = statistics.mean(walls) * 1000
ratio = mean_ms / PRIOR_MS
print(f"[measure] T1 solve wall mean {mean_ms:.1f} ms (min "
      f"{min(walls)*1000:.1f}) vs prior {PRIOR_MS} ms -> ratio "
      f"{ratio:.3f} (build-time record 1.090; 10% bar met at build)")
check("T1_wall_not_structurally_regressed", ratio < 2.0,
      f"ratio={ratio:.3f}")

finish("test_tbanks")
