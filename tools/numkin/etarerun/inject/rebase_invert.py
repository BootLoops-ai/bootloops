#!/usr/bin/env python3
"""
rebase_invert.py  (A31)

Build a rebase table  {foreign_master -> Σ c_j · M1_master_j}  purely from a
source family's kira_target.m, by inverting the within-sector basis change.

Principle:  In a sector S where a source (M89/M6/M5) chose masters {F_i} and
M1 chose {M_k} (same rank r), each M_k is in the source's target list (it's a
dotted version of some F), so the source's kira_target.m contains
      M_k = Σ_i A_{ki} F_i  +  (lower-sector source-masters)
Restricting to top-weight (sector-S) terms gives an r×r matrix A over
Q(eta) [or Q(d,eta)].  Invert:  F = A^{-1} (M − lower), then substitute
lower-sector foreigns recursively (bottom-up over sectors).

All arithmetic is exact (sympy over Q(eta)).  For a numeric-d slice, coeffs
are rational functions of eta only ⇒ small.

Output: rebase.m in kira_target.m fragment format ({F_i -> ... ,}) suitable
for extract_shared_rows.py --rebase.
"""
import os, re, sys, json, argparse
from fractions import Fraction
sys.path.insert(0, os.environ.get('ETANUMD_INJECT_DIR',
    os.path.dirname(os.path.abspath(__file__))))  # sibling-module dir (env-overridable)
from extract_shared_rows import stream_rows, parse_row, parse_int_list, fmt_int
import sympy as sp

def sec_of(idx):
    return sum((1<<i) for i in range(10) if idx[i]>0)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--m1-masters', required=True)
    ap.add_argument('--src-masters', required=True)
    ap.add_argument('--src-ktm', required=True, help="source kira_target.m")
    ap.add_argument('--out', required=True)
    ap.add_argument('--vars', default='eta',
                    help="comma-sep symbol list in coeffs (default: eta)")
    args = ap.parse_args()

    syms = {s: sp.symbols(s) for s in args.vars.split(',')}
    def coeff(s):
        return sp.nsimplify(sp.sympify(s.replace('^','**'), locals=syms),
                            rational=True)

    M1 = set(parse_int_list(args.m1_masters))
    SRC= set(parse_int_list(args.src_masters))
    foreign = sorted(SRC - M1, key=lambda i: (sec_of(i), i))
    print(f"[rebase] foreign masters (SRC\\M1): {len(foreign)}")
    fsecs = sorted(set(sec_of(f) for f in foreign))

    # Load ALL source rows keyed on lhs
    rows = {}
    for raw in stream_rows(args.src_ktm):
        pr = parse_row(raw)
        if pr is None: continue
        lhs, terms = pr
        rows[lhs] = terms
    print(f"[rebase] loaded {len(rows)} rows from src kira_target.m")

    # For each foreign sector S: collect the basis-change system.
    # F-basis = foreign masters in S; M-basis = M1 masters in S that are ALSO
    # in rows (i.e. src reduced them).
    rebase = {}   # F_idx -> list[(M1_idx, sympy_expr)]
    unresolved = []
    for S in fsecs:
        Fset = [f for f in foreign if sec_of(f)==S]
        # candidate M1-masters in this sector that src reduced
        Mset = [m for m in M1 if sec_of(m)==S and m in rows]
        # need |Mset| >= |Fset|
        if len(Mset) < len(Fset):
            print(f"[rebase] sec{S}: |Fset|={len(Fset)} > |Mset|={len(Mset)} — SKIP")
            unresolved += Fset
            continue
        # Build A (|M|×|F|) and rhs vector (M − lower):  M_k = Σ A_ki F_i + L_k
        # where L_k = Σ over (src-masters not in Fset) c·(that master)
        A = sp.zeros(len(Mset), len(Fset))
        L = [dict() for _ in Mset]   # lower[k]: {src_master_idx: sympy}
        Fpos = {f:i for i,f in enumerate(Fset)}
        for k,m in enumerate(Mset):
            for idx,c in rows[m]:
                cc = coeff(c)
                if idx in Fpos:
                    A[k,Fpos[idx]] += cc
                else:
                    L[k][idx] = L[k].get(idx, sp.Integer(0)) + cc
        # solve A^T x = e_i  for each F_i? Actually we want F = A⁺ (M − L).
        # Take first |F| rows that give invertible A (usually |M|>|F| so pick).
        # simplest: least-squares over Q(eta) is ugly; instead, pick a subset.
        # try all |F|-subsets of M until det≠0 (small: |F|≤3 typ.)
        from itertools import combinations
        solved=False
        for combo in combinations(range(len(Mset)), len(Fset)):
            Asub = A.extract(list(combo), list(range(len(Fset))))
            if Asub.det() == 0:
                continue
            Ainv = Asub.inv()
            # F_i = Σ_k Ainv[i,k] (M_{combo[k]} − L_{combo[k]})
            for i,F in enumerate(Fset):
                terms = {}
                for kk,krow in enumerate(combo):
                    a = sp.cancel(Ainv[i,kk])
                    if a == 0: continue
                    terms[Mset[krow]] = terms.get(Mset[krow], sp.Integer(0)) + a
                    for idx,cc in L[krow].items():
                        terms[idx] = terms.get(idx, sp.Integer(0)) - a*cc
                rebase[F] = terms
            solved=True
            print(f"[rebase] sec{S}: solved {len(Fset)}F via {len(Fset)}/{len(Mset)}M rows")
            break
        if not solved:
            print(f"[rebase] sec{S}: NO invertible {len(Fset)}x{len(Fset)} block — SKIP")
            unresolved += Fset

    # recursive substitution: rebase RHS may contain OTHER foreign masters
    # (from lower sectors).  bottom-up over sectors already ensures those
    # are in `rebase` — substitute.
    for F in sorted(rebase, key=lambda i: sec_of(i)):
        changed=True
        guard=0
        while changed and guard<20:
            changed=False; guard+=1
            new={}
            for idx,c in rebase[F].items():
                if idx in M1:
                    new[idx]=new.get(idx,sp.Integer(0))+c
                elif idx in rebase:
                    changed=True
                    for j,cj in rebase[idx].items():
                        new[j]=new.get(j,sp.Integer(0))+c*cj
                else:
                    # source-master ∈ M1? no. source-master not foreign? then ∈ M1.
                    # if we get here it's a foreign master not solved — leave.
                    new[idx]=new.get(idx,sp.Integer(0))+c
            rebase[F]=new
        # simplify
        rebase[F]={idx:sp.cancel(c) for idx,c in rebase[F].items() if c!=0}

    # verify: all RHS ∈ M1
    bad=0
    for F,terms in rebase.items():
        for idx in terms:
            if idx not in M1:
                bad+=1
    print(f"[rebase] built {len(rebase)} rebase rows; RHS-not-M1 residuals={bad}; "
          f"unresolved={len(unresolved)}")

    # emit as kira_target.m fragment (for extract_shared_rows --rebase)
    with open(args.out,"w") as fo:
        fo.write("{\n")
        for F,terms in rebase.items():
            fo.write(fmt_int(F)+" -> \n")
            for idx,c in terms.items():
                cs = str(sp.cancel(c)).replace('**','^')
                fo.write(f" + {fmt_int(idx)}*({cs})\n")
            fo.write(",\n")
        fo.write("}\n")
    print(f"[rebase] wrote {args.out}")
    json.dump({"n_foreign":len(foreign),"n_solved":len(rebase),
               "n_unresolved":len(unresolved),"rhs_notM1":bad},
              open(args.out+".json","w"),indent=2)

if __name__=="__main__":
    main()
