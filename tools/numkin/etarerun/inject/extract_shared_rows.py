#!/usr/bin/env python3
"""
extract_shared_rows.py  (A31)

Given one or more source kira_target.m files (η-DE reduction tables from
M5/M6/M89 slices) plus M1's master and target lists, extract the subset of
reduction rows whose LHS ∈ M1_targets and every RHS master ∈ M1_masters.

Such a row is a valid M1 η-DE row (SHARD_FEAS §1 η-coincidence rule) and can
be spliced directly into M1's kira_target.m, or fed to kira as extra_relations.

Outputs (all under --out prefix):
  <out>.rows.m         : rows in kira_target.m format ({ ... , })
  <out>.uds.kira       : rows in kira extra_relations format (blank-sep blocks)
  <out>.lhs.txt        : list of covered LHS integrals (one per line)
  <out>.stats.json     : counts per source + rejection reasons

Format of a kira_target.m row (kira2math output; <FAM> = the kira family
name, env ETANUMD_FAM):
  <FAM>[a1,...,a15] ->
   + <FAM>[b1,...,b15]*((coeff))
   ...
  ,
Rows are separated by ',', file wrapped in { ... }.  All coefficients are
exact-rational polynomials in (d,eta) or (eta) — never float.

Usage:
  extract_shared_rows.py \
      --m1-masters .../work_M1_sec1022/.../results/<FAM>/masters \
      --m1-target  .../work_M1_sec1022/.../target \
      --source M89:.../M89_eps_1_7/results/<FAM>/kira_target.m \
      [--source M6:...  --source M5:...] \
      [--diff-prop M6:4 --diff-prop M89:3]   # a_j==0 guard (SHARD_FEAS §1) \
      --out .../inject/uds_shared_1_7
"""
import os, re, sys, json, argparse
from collections import OrderedDict

FAM = os.environ.get('ETANUMD_FAM')
if not FAM:
    raise SystemExit('extract_shared_rows: set ETANUMD_FAM (kira family name)')
_F = re.escape(FAM)
INT_RE = re.compile(_F + r'\[\s*([-\d,\s]+)\s*\]')

def parse_int(s):
    m = INT_RE.search(s)
    if not m:
        return None
    return tuple(int(x) for x in m.group(1).replace(' ', '').split(','))

def parse_int_list(path):
    """read a kira integral list file (blank-line separated, # comments ok)"""
    out = []
    with open(path) as f:
        for line in f:
            line = line.split('#')[0].strip()
            if not line:
                continue
            idx = parse_int(line)
            if idx is not None:
                out.append(idx)
    return out

def fmt_int(idx):
    return FAM + "[" + ",".join(str(x) for x in idx) + "]"

# ---- kira_target.m tokenizer ------------------------------------------------
# The file is one Mathematica list { row , row , ... }.
# A row is  LHS -> ± RHS*(coef) ± RHS*(coef) ... .
# Separator between rows is a top-level ',' (paren-depth 0) that follows the
# last term of a row.  We stream-tokenize on '\n,' since kira2math always puts
# the row-separating comma on its own line (verified: M89_eps_1_7).
# For robustness we ALSO track paren depth.

def stream_rows(path):
    """yield raw row strings 'LHS -> + RHS*(c) ...' with no wrapping braces."""
    buf = []
    depth = 0
    started = False
    with open(path) as f:
        for raw in f:
            line = raw.rstrip('\n')
            s = line.strip()
            if not started:
                # skip until first '{'
                p = line.find('{')
                if p < 0:
                    continue
                line = line[p+1:]
                s = line.strip()
                started = True
            # end of file '}'
            if s == '}':
                if buf:
                    yield "".join(buf)
                return
            # paren depth of this line
            depth += line.count('(') - line.count(')')
            # a bare ',' at depth 0 terminates the row
            if s == ',' and depth == 0:
                if buf:
                    yield "".join(buf)
                    buf = []
                continue
            buf.append(line + "\n")
    if buf:
        yield "".join(buf)

TERM_RE = re.compile(
    r'([+\-])\s*' + _F + r'\[\s*([-\d,\s]+)\s*\]\s*\*\s*\((.*?)\)\s*(?=[+\-]\s*' + _F + r'\[|\Z)',
    re.DOTALL)

def parse_row(row_txt):
    """return (lhs_idx, [(rhs_idx, coeff_str_signed), ...]) or None"""
    # split on '->'
    p = row_txt.find('->')
    if p < 0:
        return None
    lhs = parse_int(row_txt[:p])
    rhs_txt = row_txt[p+2:].strip()
    if lhs is None:
        return None
    # kira2math emits ' + I*(c)\n + I*(c)...' — leading sign always present.
    # Use TERM_RE with a lookahead terminator so nested parens in coeff are
    # captured greedily up to the next ' + <FAM>[' or EOS.  This is safe
    # because kira never emits the family token inside a coefficient.
    terms = []
    for m in TERM_RE.finditer(rhs_txt):
        sign, idx_s, coeff = m.group(1), m.group(2), m.group(3).strip()
        idx = tuple(int(x) for x in idx_s.replace(' ', '').split(','))
        # bake sign into coeff so uds emit is trivially correct
        if sign == '-':
            coeff = f"-({coeff})"
        terms.append((idx, coeff))
    if not terms:
        # a row with no RHS terms means LHS is a master (kira sometimes emits
        # 'I -> 0' or just 'I' — treat as LHS=0 or skip; here we skip: no info)
        return (lhs, [])
    return (lhs, terms)

# ---- main -------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--m1-masters', required=True)
    ap.add_argument('--m1-target',  required=True)
    ap.add_argument('--source', action='append', required=True,
                    help="TAG:path/to/kira_target.m (repeat)")
    ap.add_argument('--diff-prop', action='append', default=[],
                    help="TAG:j — require a_j==0 for LHS+RHS from that source")
    ap.add_argument('--out', required=True)
    ap.add_argument('--reuse-map', default=None,
                    help="M1_reuse_map.json (informational cross-check only)")
    ap.add_argument('--rebase', default=None,
                    help="kira_target.m expressing foreign masters onto M1 masters; "
                         "substituted into RHS before the M1_masters check")
    ap.add_argument('--emit-foreign', default=None,
                    help="write list of RHS integrals ∉ M1_masters (for building --rebase)")
    args = ap.parse_args()

    M1_masters = set(parse_int_list(args.m1_masters))
    M1_target  = set(parse_int_list(args.m1_target))

    # rebase table: foreign_idx -> [(m1_idx, coeff), ...]
    rebase = {}
    if args.rebase:
        for raw in stream_rows(args.rebase):
            pr = parse_row(raw)
            if pr is None:
                continue
            lhs, terms = pr
            if all(t[0] in M1_masters for t in terms):
                rebase[lhs] = terms
        print(f"[extract] rebase table: {len(rebase)} foreign→M1 rows")
    diffprop = {}
    for dp in args.diff_prop:
        tag, j = dp.split(':')
        diffprop.setdefault(tag, []).append(int(j))

    print(f"[extract] M1 masters={len(M1_masters)}  M1 targets={len(M1_target)}")

    # OrderedDict keyed on lhs_idx to dedupe across sources (first wins)
    kept = OrderedDict()
    foreign = set()   # RHS integrals ∉ M1_masters seen anywhere
    stats = {"sources": {}, "total_kept": 0, "M1_targets": len(M1_target),
             "M1_masters": len(M1_masters), "n_rebase": len(rebase)}

    for src in args.source:
        tag, path = src.split(':', 1)
        st = {"rows_in_file": 0, "lhs_in_M1target": 0,
              "rhs_all_M1masters": 0, "diffprop_violations": 0,
              "kept_new": 0, "kept_dup": 0, "empty_rhs": 0, "rebased": 0}
        dprops = diffprop.get(tag, [])
        for raw in stream_rows(path):
            pr = parse_row(raw)
            if pr is None:
                continue
            lhs, terms = pr
            st["rows_in_file"] += 1
            if not terms:
                st["empty_rhs"] += 1
                continue
            # diff-prop guard: for M6/M89, LHS and every RHS must have a_j==0
            if dprops:
                if any(lhs[j-1] != 0 for j in dprops) or \
                   any(t[0][j-1] != 0 for t in terms for j in dprops):
                    st["diffprop_violations"] += 1
                    continue
            if lhs not in M1_target:
                continue
            st["lhs_in_M1target"] += 1
            # rebase step: substitute foreign RHS via rebase table
            new_terms = []
            did_rebase = False
            ok = True
            for idx, coeff in terms:
                if idx in M1_masters:
                    new_terms.append((idx, coeff))
                elif idx in rebase:
                    did_rebase = True
                    for ridx, rcoeff in rebase[idx]:
                        new_terms.append((ridx, f"({coeff})*({rcoeff})"))
                else:
                    foreign.add(idx)
                    ok = False
            if not ok:
                continue
            st["rhs_all_M1masters"] += 1
            if did_rebase:
                st["rebased"] += 1
            # collect like-terms (multiple RHS may map to same M1-master after rebase)
            from collections import defaultdict as dd
            acc = dd(list)
            for idx, coeff in new_terms:
                acc[idx].append(coeff)
            terms_final = [(idx, coeffs[0] if len(coeffs) == 1
                                 else "+".join(f"({c})" for c in coeffs))
                           for idx, coeffs in acc.items()]
            if lhs in kept:
                st["kept_dup"] += 1
                continue
            kept[lhs] = (tag, terms_final)
            st["kept_new"] += 1
        stats["sources"][tag] = st
        print(f"[extract] {tag}: {st}")

    stats["total_kept"] = len(kept)
    stats["M1_targets_residual"] = len(M1_target) - len(kept)

    # ---- emit .rows.m (kira_target.m fragment) -----------------------------
    with open(args.out + ".rows.m", "w") as f:
        for lhs, (tag, terms) in kept.items():
            f.write(fmt_int(lhs) + " -> \n")
            for idx, coeff in terms:
                f.write(f" + {fmt_int(idx)}*({coeff})\n")
            f.write(",\n")

    # ---- emit .uds.kira (extra_relations format: sum=0 blocks) -------------
    # equation:  LHS - Σ c_j RHS_j = 0
    with open(args.out + ".uds.kira", "w") as f:
        for lhs, (tag, terms) in kept.items():
            f.write(f"{fmt_int(lhs)}*(-1)\n")
            for idx, coeff in terms:
                f.write(f"{fmt_int(idx)}*({coeff})\n")
            f.write("\n")

    # ---- emit .lhs.txt -----------------------------------------------------
    with open(args.out + ".lhs.txt", "w") as f:
        for lhs in kept:
            f.write(fmt_int(lhs) + "\n")

    # ---- residual target file (M1 target minus covered LHS) ----------------
    with open(args.out + ".target_residual", "w") as f:
        for t in parse_int_list(args.m1_target):  # preserve order
            if t not in kept:
                f.write(FAM + "[" + ", ".join(str(x) for x in t) + "]\n\n")

    stats["n_foreign_unresolved"] = len(foreign)
    with open(args.out + ".stats.json", "w") as f:
        json.dump(stats, f, indent=2)

    if args.emit_foreign:
        with open(args.emit_foreign, "w") as f:
            for idx in sorted(foreign):
                f.write(FAM + "[" + ", ".join(str(x) for x in idx) + "]\n\n")
        print(f"[extract] wrote {len(foreign)} unresolved foreign masters -> {args.emit_foreign}")

    print(f"[extract] KEPT {len(kept)} shared rows;  "
          f"M1 residual targets = {stats['M1_targets_residual']}")
    print(f"[extract] wrote: {args.out}.{{rows.m,uds.kira,lhs.txt,target_residual,stats.json}}")

if __name__ == "__main__":
    main()
