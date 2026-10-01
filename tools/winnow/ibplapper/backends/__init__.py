"""backends — optional accelerated eliminator backends.

The default ibplapper path (backend='cpu') is the verbatim-lifted sparse
dict engine (engine.py) and is UNTOUCHED by this package: importing
`ibplapper` never imports a backend (torch stays a lazy, optional
dependency; api.eliminate imports backends.stratified only when
backend='dense' is requested).

Modules:
  dense       eliminate_dense — dense int64 torch eliminator, contract-
              identical to engine.eliminate_fast on a stratum block
              (device='cpu' fully supported; device='cuda' same code path).
  stratified  stratified_solve_dense — mirror of engine.stratified_solve
              (B2F/B2FT) dispatching each stratum block to eliminate_dense
              under a dense-buffer budget, else engine.eliminate_fast.

Import explicitly:  from ibplapper.backends.dense import eliminate_dense
(no submodule is imported here, so `import ibplapper.backends` stays free
of the torch import).
"""
