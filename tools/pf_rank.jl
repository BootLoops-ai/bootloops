# ============================================================================
# pf_rank.jl — COMPAT SHIM, FOREVER.
# Canonical member: tools/dipstick/pf_rank.jl (verb `order` of DIPSTICK, the
# unified pre-compute triage package; manual: tools/dipstick/GUIDE.md).
# include() of THIS path keeps working for the existing call form —
# maxcut/decisive.py Track B (PF_RANK_JL = tools/pf_rank.jl, `julia -e
# include(...)`) and any caller that includes the flat path directly —
# defining pf_probe / pf_probe_gen / pf_probe_esc AND the internals (_RE,
# _feed!, _mm, _row, _pf_core) in the caller's module, exactly as before.
# (The validation harness tools/dipstick/pf_rank_validate.jl includes the
# dipstick member directly, not this shim.)
# Do not add logic here; see tools/dipstick/pf_rank.jl.
# ============================================================================
include(joinpath(@__DIR__, "dipstick", "pf_rank.jl"))
