#!/usr/bin/env python3
"""kira_parse.py — pure-regex parser for Kira's kira2math output.

Self-contained on purpose: no module-level file loads, no imports beyond re.
Returns { lhs_tuple : [ (rhs_tuple, coef_string), ... ] }; coef_string is the
raw rational function in (d, var[, s, t]).
"""
import re

_FAMINT = re.compile(r'([A-Za-z_][A-Za-z0-9_]*)\[\s*([0-9,\-\s]+)\]')


def _idx(tok):
    return tuple(int(x) for x in tok.replace(' ', '').split(','))


def parse_kira2math(path, family=None):
    txt = open(path).read().strip()
    if txt.startswith('{'):
        txt = txt[1:]
    if txt.endswith('}'):
        txt = txt[:-1]
    # kira2math emits the inter-rule comma on its own line
    chunks, cur = [], []
    for ln in txt.splitlines():
        if ln.strip() == ',':
            chunks.append('\n'.join(cur)); cur = []
        else:
            cur.append(ln)
    if cur:
        chunks.append('\n'.join(cur))
    out = {}
    for ch in chunks:
        ch = ch.strip()
        if not ch:
            continue
        m = _FAMINT.match(ch)
        if not m:
            continue
        fam, idx_s = m.groups()
        if family is not None and fam.lower() != family.lower():
            continue
        lhs = _idx(idx_s)
        rest = ch[m.end():].lstrip()
        if not rest.startswith('->'):
            continue          # master mapped to itself (no arrow)
        rest = rest[2:]
        terms = []
        pos = 0; n = len(rest)
        while pos < n:
            while pos < n and rest[pos] in ' \n\r\t+':
                pos += 1
            if pos >= n:
                break
            mm = _FAMINT.match(rest, pos)
            if not mm:
                break
            rhs = _idx(mm.group(2)); pos = mm.end()
            while pos < n and rest[pos] in ' \n\r\t':
                pos += 1
            if pos < n and rest[pos] == '*':
                pos += 1
            while pos < n and rest[pos] in ' \n\r\t':
                pos += 1
            if pos < n and rest[pos] == '(':
                depth = 0; start = pos
                while pos < n:
                    c = rest[pos]
                    if c == '(':
                        depth += 1
                    elif c == ')':
                        depth -= 1
                        if depth == 0:
                            pos += 1; break
                    pos += 1
                coef = rest[start:pos]
            else:                      # bare coefficient (rare)
                nm = _FAMINT.search(rest, pos)
                end = nm.start() if nm else n
                coef = rest[pos:end].strip().rstrip('+').strip()
                pos = end
            terms.append((rhs, coef))
        out[lhs] = terms
    return out
