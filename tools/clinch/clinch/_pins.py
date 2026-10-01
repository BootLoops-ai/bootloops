"""clinch._pins — pinned source sha256-16 of every clinch module.

Regenerate ONLY as part of a deliberate tool change (manual updated and
battery green in the same change):
  python3 -c "import clinch._pins as p; p._regen()"
"""
PINS = {
    "jets.py": "ab1e53d151d5af9e",
    "oracle_v36.py": "ef21775eb5e13c39",
    "oracle_v31.py": "64e0baff67b67018",
    "adapter_v31.py": "50bbd33c8be6c6f6",
    "engine.py": "905b5c699767a014",
    "kkt_pd.py": "af5b56bc7aa08c3f",
    "cert.py": "c186f18b1b04eef4",
}


def _regen():
    import hashlib, os, re
    here = os.path.dirname(os.path.abspath(__file__))
    src = open(os.path.join(here, "_pins.py")).read()
    for m in list(PINS):
        h = hashlib.sha256(open(os.path.join(here, m), "rb").read()).hexdigest()[:16]
        src = re.sub(r'"%s": "[0-9a-f]{16}"' % re.escape(m), '"%s": "%s"' % (m, h), src)
    open(os.path.join(here, "_pins.py"), "w").write(src)
    print("pins regenerated")
