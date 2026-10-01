# amflow_smoke fixtures
Job specs and reference output for the amflow-cpp smoke test.
- in_vac2Bprobe.json : 2L vacuum probe, 5 dotted ints.
- in_m3soft.json     : 1-loop triangle job spec, first 2 integrals (pre-cut; smoke's [:2] is idempotent).
- out_vac2Bprobe.json: vacuum reference output; agrees to full printed precision (>120 digits) across independent builds.
