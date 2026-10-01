"""CLI: python -m mixalot.cli counts.json [--gmax 3]

counts.json SCHEMA (exactly one of):
  [3, 1, 4, 1, 5]                      a bare count vector (integers), or
  {"U": [3,1,4,1,5], "gmax": 3}        object form; "U" required, "gmax" optional.
Counts must be integers (non-integer counts are REJECTED, never truncated).
Prints the verdict block (exact values as strings) + the advisory."""
import json
import sys
from fractions import Fraction

from .advisory import advise
from .core import gstar


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    gmax = 3
    if "--gmax" in argv:
        i = argv.index("--gmax")
        gmax = int(argv[i + 1])
        del argv[i:i + 2]
    if not argv:
        print(__doc__)
        return 2
    spec = json.load(open(argv[0]))
    if isinstance(spec, dict):
        if "U" not in spec:
            print("ERROR: object-form counts.json requires the key \"U\" "
                  "(see --help schema)", file=sys.stderr)
            return 2
        U = spec["U"]
        if "gmax" in spec:
            gmax = int(spec["gmax"])
    elif isinstance(spec, list):
        U = spec
    else:
        print("ERROR: counts.json must be a list or an object with \"U\"",
              file=sys.stderr)
        return 2
    r = gstar(U, gmax=gmax)
    block = {
        "U": list(map(int, U)), "k": len(U), "N": int(sum(U)), "gmax": gmax,
        "posterior": {str(g): str(p) + f"  (~{float(p):.6g})"
                      for g, p in r["posterior"].items()},
        "map_g": r["map_g"],
        "p_ge2": str(r["p_ge2"]) + f"  (~{float(r['p_ge2']):.6g})",
        "evidence_Z": {str(g): str(z) for g, z in r["evidence"].items()},
        "register": ("EXACT rational posterior under the vendored gated "
                     "evaluator (Dir(1) convention, uniform prior over "
                     f"g=1..{gmax} unless supplied); fences in README"),
        "advisory": advise(N=int(sum(U))),
    }
    print(json.dumps(block, indent=1, default=str))
    return 0


if __name__ == "__main__":
    sys.exit(main())
