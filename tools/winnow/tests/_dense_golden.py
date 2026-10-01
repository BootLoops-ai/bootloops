"""Shared golden T-bank byte-equal gates for the dense backend — the
single implementation used by tests/test_dense_backend.py (device cpu, in
the battery); a cross-device determinism harness can run the same gates
with --device cuda. Import, never duplicate.

Gates (identical names in both callers):
  GOLDEN_k2disp_T_byte_equal::dense_B2F
  GOLDEN_k2disp_T_byte_equal::dense_B2FT
  GOLDEN_lp1disp_T_byte_equal::dense_B2FT
Each replays all 4 archived slices through the ibplapper API with
backend='dense' on the requested device and compares the produced T .bin
sha256 to the archived golden bank.
"""
import hashlib
import os

import ibplapper as lap
from ibplapper.adapters import kira
from _common import D4OUT, NOMAT, SLICES4, STAGEPACKS


def sha(path):
    return hashlib.sha256(open(path, "rb").read()).hexdigest()


def replay_dense(gen, fam, policy, base, device="cpu"):
    """Replay the 4 archived slices through the API (backend='dense') into a
    fresh Bank at `base`; asserts zero leftover and all-dense stratum paths
    (a budget fallback would not witness the dense kernel)."""
    for ext in (".bin", ".idx.json"):
        if os.path.exists(base + ext):
            os.remove(base + ext)
    for si, (p, d0, e0) in enumerate(SLICES4):
        system, _sysd = kira.load(gen, fam, p, d0, e0)
        bank = lap.Bank(base, meta=None if si else {"family": fam})
        res = lap.eliminate(system, lap.Schedule(
            policy=policy, bank=bank, bank_cell="tiny",
            bank_node=f"p{p}_d{d0}_e{e0}"), backend="dense", device=device)
        assert not res.leftover
        assert all(s["path"] == "dense" for s in res.stats_per.values()), \
            "budget fallback fired on a golden fixture — gate not dense"
    return base + ".bin"


def golden_gates(out_dir, device="cpu"):
    """Yield (gate_name, passed) for the three golden byte-equal gates."""
    gen_k2 = os.path.join(NOMAT, "gen_lbl3m2L_k2disp")
    banked_k2 = os.path.join(D4OUT, "lbl3m2L_k2disp", "T_lbl3m2L_k2disp.bin")
    for pol in ("B2F", "B2FT"):
        got = replay_dense(gen_k2, "lbl3m2L_k2disp", pol,
                           os.path.join(out_dir, f"T_k2disp_dense_{pol}"),
                           device=device)
        yield f"GOLDEN_k2disp_T_byte_equal::dense_{pol}", \
            sha(got) == sha(banked_k2)
    banked_lp = os.path.join(D4OUT, "lbl3m2L_lp1disp",
                             "T_lbl3m2L_lp1disp.bin")
    got = replay_dense(os.path.join(STAGEPACKS, "gen_lbl3m2L_lp1disp"),
                       "lbl3m2L_lp1disp", "B2FT",
                       os.path.join(out_dir, "T_lp1disp_dense"),
                       device=device)
    yield "GOLDEN_lp1disp_T_byte_equal::dense_B2FT", sha(got) == sha(banked_lp)
