# repertoire_g1 — provenance

Code-only absorption of a Gate-1 complete-census certificate engine from
an external immune-repertoire study (De Boer-Perelson clone-niche
competition family): g1_cert.py (certified layer), g1_census.py
(per-cell census driver), g1_model.py (model family + anchored grid),
g1_run.py (bench/batch/assemble driver). PROVENANCE.sha256 lists the
sha256 of every shipped file.

NOT shipped: the study's receipts (out/ cells, GATE1_VERDICT.json, done
markers). The independent-arithmetic verifier (a non-arb fixed-point
interval twin used for rederivation) is declared under baller, not here.

Receipts writing: g1_run.py writes to $G1_OUT, defaulting to ./g1_out
under the caller's cwd — never inside the package tree. Set G1_OUT
explicitly for batch runs.
