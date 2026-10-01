#!/usr/bin/env python3
"""paper_seams.py — comment-swallow seam gate for a LaTeX paper tree (emitall's document side).

usage: python3 tools/emitall/paper_seams.py <master.tex> [--snapshot <dir>] [--marker <label> ...]
                                            [--allow <phrase> ...] [--quiet]
  --allow <phrase>: a legitimate lowercase sentence head (e.g. a tool name 'txburst, the backend') --
  suppresses S2 fragments containing it AND S1/S1b/S3 strikes on a live line that begins with it.
  --marker <label>: a label that identifies a fold comment block for S3 (repeatable; default 'FOLD').
exit 0 = clean; exit 1 = at least one tex-confirmed swallow; exit 2 = usage (no master, -h/--help).

A comment swallow: an edit drops a % comment block INSIDE a sentence, so the head of the sentence is
stranded on the comment's last line and the rendered text reads "...end. lowercase fragment...". The
cure is a MOVE of the comment to its own line boundary, never a rewording of the prose.

Four rows, all over the master's \\input/\\include roster (the rendered surface only):
  S1  a live line that ENDS a sentence (. ? !), then one or more % comment lines, then a live line
      that BEGINS lowercase -- the head of the sentence is stranded on the comment's last line.
  S1b a comment block whose LAST line ends on a capitalized sentence-head word with no terminal
      punctuation, followed by a lowercase live line (a shape with no legitimate counterpart).
  S2  rendered 'word. lowercase' fragments (pdftotext -layout beside the master; numeric enders
      included; abbreviations excluded; bibliography excluded) that the de-commented tex ALSO carries
      adjacently -- reading-order artifacts (dropped section glyphs, ligatures, table layout) are INFO.
  S3  (with --snapshot) a fold comment block (one carrying a --marker label) that ends on a bare
      lowercase word, or carries after its own sentence words that were LIVE text in the pre-fold
      snapshot, when the next live line begins lowercase.
  S1c (with --snapshot) every clause (>= 4 words, split at . ; : ( )) of every LIVE snapshot line
      that the edit deleted or rewrote must survive somewhere in the current file -- as live text or
      inside a comment (superseded text is conventionally kept, quoted, in a comment). A clause that
      vanished from both is a swallow no other row can see: e.g. a physical line comment-blanked to
      preserve a superseded clause that ALSO carried live words of the next sentence. A clause cut
      mid-way by an edit boundary passes if it survives as up to three contiguous word-runs.

Tabular rows (' & ') and macro lines (leading backslash) legitimately begin lowercase and are skipped.
Without a PDF beside the master (or without pdftotext) S2 is skipped and says so; S1/S1b/S3/S1c run on
the tex alone (S3 and S1c only with --snapshot). Run it on every paper master before the built document goes out, beside battery.py's
claims rows; the fixture battery lives in tests/test_paper_seams.py over specs/seams/.
"""
import glob, os, re, shutil, subprocess, sys

ABBREV = {'e.g', 'i.e', 'cf', 'vs', 'al', 'eq', 'eqs', 'fig', 'figs', 'sec', 'secs', 'ref', 'refs',
          'no', 'vol', 'pp', 'ch', 'approx', 'resp', 'etc', 'viz', 'ibid', 'op', 'loc', 'ca', 'wrt'}
UNIT_NEXT = {'and', 'or', 'to', 'of', 'in', 'at', 'by', 'vs', 'per', 'mbar', 'ma', 'ga', 'kyr', 'myr',
             'gyr', 'pp', 'cm', 'mm', 'km', 'nm', 'mol', 'au', 'd', 's', 'h', 'min', 'digits'}
MARKERS = ('FOLD',)   # default fold-block labels for S3; extend with --marker <label>


def roster(master):
    base = os.path.dirname(os.path.abspath(master))
    try:
        txt = open(master, encoding='utf-8').read()
    except OSError:
        return [master]
    live = '\n'.join(l for l in txt.split('\n') if not l.strip().startswith('%'))
    files = [master]
    for r in re.findall(r'\\(?:input|include)\{([^}]+)\}', live):
        cand = os.path.join(base, r if r.endswith('.tex') else r + '.tex')
        if os.path.exists(cand) and cand not in files:
            files.append(cand)
            # one level of nesting (a section file that \inputs sub-files)
            try:
                sub = open(cand, encoding='utf-8').read()
            except OSError:
                continue
            sub_live = '\n'.join(l for l in sub.split('\n') if not l.strip().startswith('%'))
            for r2 in re.findall(r'\\(?:input|include)\{([^}]+)\}', sub_live):
                c2 = os.path.join(base, r2 if r2.endswith('.tex') else r2 + '.tex')
                if os.path.exists(c2) and c2 not in files:
                    files.append(c2)
    return files


def s1_s1b(files, base_dir, allow=()):
    fails = []
    for path in files:
        try:
            lines = open(path, encoding='utf-8').read().split('\n')
        except OSError:
            continue
        rel = os.path.relpath(path, base_dir)
        last_live, last_comment, in_comment = '', '', False
        for i, raw in enumerate(lines, 1):
            s = raw.strip()
            if not s:
                continue
            if s.startswith('%'):
                in_comment = True
                last_comment = raw
                continue
            if in_comment:
                in_comment = False
                if any(s.startswith(a) for a in allow):
                    # an allowed lowercase sentence head (a tool name such as 'txburst, the backend')
                    last_live = raw
                    continue
                if ' & ' in s or s.startswith('\\'):
                    # a tabular row / a macro line legitimately begins lowercase
                    last_live = raw
                    continue
                prev = re.sub(r'\s*%.*$', '', last_live).rstrip().rstrip('}').rstrip()
                if prev and prev[-1] in '.?!' and re.match(r'[a-z]', s):
                    fails.append(f'S1 {rel}:{i} live line begins lowercase after a comment block that follows a sentence end: "{s[:60]}"')
                lc = last_comment.strip().lstrip('%').strip()
                toks = lc.split()
                if (toks and re.match(r'[a-z]', s) and not re.search(r'[.:;?!)\]"\'`}]$', lc)
                        and re.match(r'^[A-Z][a-z]+$', toks[-1])):
                    fails.append(f'S1b {rel}:{i} comment ends with the sentence head "{" ".join(toks[-3:])}" and the live line begins lowercase: "{s[:50]}"')
            last_live = raw
    return fails


def s2(master, files, allow):
    pdf = os.path.splitext(master)[0] + '.pdf'
    if not os.path.exists(pdf) or not shutil.which('pdftotext'):
        return [], 'S2 skipped (no PDF or no pdftotext)'
    try:
        txt = subprocess.run(['pdftotext', '-layout', '-enc', 'UTF-8', pdf, '-'],
                             capture_output=True, text=True, timeout=180).stdout
    except (subprocess.SubprocessError, OSError) as e:
        return [f'S2 pdftotext failed ({e})'], 'S2 failed'
    m = re.search(r'\n\s*References\s*\n', txt)
    body = txt[:m.start()] if m else txt
    flat = re.sub(r'\s+', ' ', body)
    tex_flat = []
    for path in files:
        try:
            raw = open(path, encoding='utf-8').read()
        except OSError:
            continue
        live = '\n'.join(l for l in raw.split('\n') if not l.strip().startswith('%'))
        live = re.sub(r'(?<!\\)%.*', '', live).replace('~', ' ')
        tex_flat.append(re.sub(r'\s+', ' ', live))
    tex_flat = ' '.join(tex_flat)
    fails, hits, artifacts = [], 0, 0
    for mm in re.finditer(r'([A-Za-z][a-z]{1,}|\d+)\. ([a-z]{2,})', flat):
        word, nxt = mm.group(1), mm.group(2)
        if word.lower() in ABBREV or (word.isdigit() and nxt in UNIT_NEXT):
            continue
        frag = flat[max(0, mm.start() - 30):mm.end() + 30]
        hits += 1
        if any(a in frag for a in allow):
            continue
        if not re.search(re.escape(word) + r'\.\s+' + re.escape(nxt) + r'\b', tex_flat):
            artifacts += 1
            continue
        fails.append(f'S2 rendered lowercase-after-period confirmed in tex: "...{frag}..."')
    return fails, f'S2: {hits} candidates, {artifacts} rendering artifacts (INFO), {len(fails)} tex-confirmed'


def s3(files, base_dir, snapshot_dir, allow=(), markers=MARKERS):
    fails = []
    snap_live = ''
    if snapshot_dir and os.path.isdir(snapshot_dir):
        for f in files:
            rel = os.path.relpath(f, base_dir)
            for cand in (os.path.join(snapshot_dir, rel), os.path.join(snapshot_dir, os.path.basename(f))):
                if os.path.exists(cand):
                    t = open(cand, encoding='utf-8').read()
                    snap_live += ' ' + re.sub(r'\s+', ' ', '\n'.join(l for l in t.split('\n') if not l.strip().startswith('%')))
                    break
    for path in files:
        lines = open(path, encoding='utf-8').read().split('\n')
        rel = os.path.relpath(path, base_dir)
        block, start = [], None
        for i, raw in enumerate(lines + [''], 1):
            st = raw.strip()
            if st.startswith('%'):
                if start is None:
                    start = i
                block.append(st)
                continue
            if block:
                joined = ' '.join(block)
                last = block[-1].lstrip('%').strip()
                next_lower = bool(re.match(r'[a-z]', st)) and ' & ' not in st and not any(st.startswith(a) for a in allow)
                if next_lower and any(mk in joined for mk in markers) and re.search(r'[a-z]$', last) and not re.search(r'[.:;)\]}\'"`]\s*$', last):
                    fails.append(f'S3 {rel}:{i-1} fold comment block ends on a bare word ("...{last[-50:]}") before a lowercase live line')
                if snap_live and next_lower and any(mk in joined for mk in markers):
                    for k, cl in enumerate(block):
                        bodyl = re.sub(r'^%+', '', cl).strip()
                        if bodyl.count('"') % 2 == 1:
                            continue
                        m = re.search(r'[.)]\s+([a-z][^.]{3,}?)$', bodyl)
                        if m and len(m.group(1).split()) >= 2 and m.group(1).strip() in snap_live:
                            fails.append(f'S3 {rel}:{start+k} comment tail was LIVE text in the pre-fold snapshot ("...{m.group(1)[-60:]}")')
                block, start = [], None
    return fails


def _s1c_pieces(c, cur_all, maxpieces=3, minw=2):
    """S1c second chance: a snapshot clause that an edit boundary cuts mid-clause survives as up to
    `maxpieces` contiguous word-runs (head live / middle quoted / tail live), each >= `minw` words and
    each present in the live-or-quoted stream. True = covered."""
    w = c.split(); n = len(w)
    have = lambda i, j: j - i >= minw and ' '.join(w[i:j]) in cur_all
    for i in range(minw, n - minw + 1):
        if have(0, i) and have(i, n):
            return True
    if maxpieces >= 3:
        for i in range(minw, n):
            if not have(0, i):
                continue
            for j in range(i + minw, n - minw + 1):
                if have(i, j) and have(j, n):
                    return True
    return False


def s1c(files, base_dir, snapshot_dir):
    """S1c: clauses of deleted/rewritten LIVE snapshot lines that survive neither as live text nor inside
    a comment of the current file. Needs --snapshot = the pre-edit copy of the tree."""
    import difflib
    fails = []
    if not (snapshot_dir and os.path.isdir(snapshot_dir)):
        return fails
    norm = lambda x: ' '.join(x.split())
    for path in files:
        rel = os.path.relpath(path, base_dir)
        snap_path = None
        for cand in (os.path.join(snapshot_dir, rel), os.path.join(snapshot_dir, os.path.basename(path))):
            if os.path.exists(cand):
                snap_path = cand
                break
        if not snap_path:
            continue
        try:
            snap = open(snap_path, encoding='utf-8').read().split('\n')
            cur_raw = open(path, encoding='utf-8').read()
        except OSError:
            continue
        cur_lines = cur_raw.split('\n')
        # everything the current file carries, live AND comment text, whitespace-normalized; macro
        # openers (\caption{, \emph{ ...) and braces are dropped on both sides so a clause that begins
        # inside a macro argument still matches its quoted copy
        demac = lambda x: re.sub(r'[{}]', ' ', re.sub(r'\\[A-Za-z]+\*?\{', ' ', x))
        cur_all = norm(demac(' '.join(re.sub(r'^\s*%+\s?', '', l) for l in cur_lines)))
        sm = difflib.SequenceMatcher(None, snap, cur_lines, autojunk=False)
        for tag, i1, i2, _j1, _j2 in sm.get_opcodes():
            if tag not in ('delete', 'replace'):
                continue
            for k in range(i1, i2):
                line = snap[k]
                if not line.strip() or line.strip().startswith('%'):
                    continue
                live = re.sub(r'(?<!\\)%.*$', '', line)
                for clause in re.split(r'[.;:()]', demac(live)):
                    c = norm(clause)
                    if len(c.split()) < 4:
                        continue
                    if c not in cur_all and not _s1c_pieces(c, cur_all):
                        fails.append(f'S1c {rel}: snapshot line {k+1} live clause vanished in the edit (neither live nor quoted): "{c[:80]}"')
    return fails


def main(argv):
    if len(argv) < 2 or argv[1] in ('-h', '--help'):
        print(__doc__); return 2
    master = os.path.abspath(argv[1])
    snap = None; allow = []; quiet = False; markers = []
    i = 2
    while i < len(argv):
        if argv[i] == '--snapshot': snap = argv[i + 1]; i += 2
        elif argv[i] == '--allow': allow.append(argv[i + 1]); i += 2
        elif argv[i] == '--marker': markers.append(argv[i + 1]); i += 2
        elif argv[i] == '--quiet': quiet = True; i += 1
        else: i += 1
    base_dir = os.path.dirname(master)
    files = roster(master)
    fails = s1_s1b(files, base_dir, allow)
    f2, note = s2(master, files, allow)
    fails += f2
    fails += s3(files, base_dir, snap, allow, tuple(markers) or MARKERS)
    fails += s1c(files, base_dir, snap)
    if not quiet:
        for f in fails: print('FAIL', f)
        print(f'INFO {note}; roster {len(files)} files')
    print(f'paper_seams: {"CLEAN" if not fails else str(len(fails)) + " swallow(s)"} -- {os.path.basename(master)}')
    return 1 if fails else 0


if __name__ == '__main__':
    sys.exit(main(sys.argv))
