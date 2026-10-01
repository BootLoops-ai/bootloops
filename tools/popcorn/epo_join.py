#!/usr/bin/env python3
"""epo_join — ancestral-state join of a sites table against an EPO
ancestral-allele FASTA (polarization by outgroup-reconstructed ancestor).

WHAT IT COMPUTES
  Streams a table of variant sites (chrom, pos, ref, alt[, ...] or
  pos, ref, alt[, ...]) against a per-chromosome ancestral-allele FASTA of
  the Ensembl EPO kind and, for every site, looks up the ancestral base
  anc = FASTA[pos] (1-based; '.' when pos lies outside the sequence),
  classifies the call by confidence, and writes the polarized table

      pos  ref  alt  anc  conf          (tab-separated, one header row)

  with conf in {'high', 'low', 'none'}. Alongside the table it accumulates
  integer counters: total sites, sites outside the FASTA, sites of other
  chromosomes skipped, the confidence histogram over all sites, and, over
  biallelic SNVs only (len(ref) == len(alt) == 1, both in ACGT), the
  orientation counters
      high/low x {anc_eq_ref, anc_eq_alt, anc_neither},  unpolarizable
  i.e. whether the ancestral base equals REF (REF ancestral, ALT derived),
  equals ALT (ALT ancestral, REF derived), or neither (a mismatch class that
  cannot be polarized by this ancestor), and how many SNVs carry no call.
  polarization_rates() turns those counters into the usual summary
  fractions (mismatch rate, low-confidence rate, unpolarizable fraction).
  Counting only: nothing here fits a model.

CONVENTION (Ensembl EPO ancestral-allele FASTA; the same case convention is
carried by the AA field of the 1000 Genomes Project releases). Ortheus
infers the ancestral state from the EPO alignment and grades it by comparing
the call with the sister sequence and with the next-deeper ancestor:
  A C G T   high-confidence call: ancestral state supported by the other
            two sequences
  a c g t   low-confidence call: ancestral state supported by one sequence
            only
  N         failure: ancestral state not supported by any other sequence
  -         the extant species carries an insertion at this position
  .         no alignment coverage
  classify(anc) maps these to 'high' | 'low' | 'none' ('N', 'n', '-', '.',
  and anything else, are 'none'). Only the FIRST record of the FASTA is read
  (one chromosome per file, as Ensembl ships them); further records are
  ignored and flagged in the counters.

INPUT LAYOUTS (sites table; '#'-prefixed lines are skipped everywhere;
'.gz' accepted; '-' reads standard input)
  chrom    columns chrom, pos, ref, alt[, ...]  (e.g. `cut -f1,2,4,5` of a
           VCF body; a header row whose second name is 'pos', any case, is
           consumed). With chrom="chr21" rows whose first column is neither
           'chr21' nor '21' are skipped and counted (skipped_other_chrom);
           with chrom=None every row is used.
  nochrom  columns pos, ref, alt[, ...] with an optional header row whose
           first three names are exactly  pos ref alt.
  auto     (default) picks 'nochrom' when the first data line starts with
           the header token 'pos' or has an integer in column 1 and fewer
           than four columns or a non-integer column 2; otherwise 'chrom'.

REGISTER: DATA-UTILITY. Integer counters and string columns; the only floats
are the summary fractions, plain x/d of the integer counters. Nothing here
is exact-algebraic, certified, or an estimator; do not quote it beside the
certified registers of the package as if it were one. Comparisons of the
output table across runs should be made on the decompressed payload (a
'.gz' output carries an mtime in its header).

REFERENCES
  Paten, Herrero, Beal, Fitzgerald, Birney (2008), "Enredo and Pecan:
    genome-wide mammalian consistency-based multiple alignment with
    paralogs", Genome Res. 18:1814-1828, doi:10.1101/gr.076554.108.
  Paten, Herrero, Fitzgerald, Beal, Flicek, Holmes, Birney (2008),
    "Genome-wide nucleotide-level mammalian ancestor reconstruction",
    Genome Res. 18:1829-1843, doi:10.1101/gr.076521.108  (Ortheus; the
    EPO = Enredo-Pecan-Ortheus ancestral calls).
  Herrero et al. (2016), "Ensembl comparative genomics resources",
    Database 2016:bav096, doi:10.1093/database/bav096.
  The 1000 Genomes Project Consortium (2015), "A global reference for human
    genetic variation", Nature 526:68-74, doi:10.1038/nature15393
    (ancestral-allele annotation with the upper/lower-case convention).

USAGE
  from epo_join import join_sites, polarization_rates
  counts = join_sites("homo_sapiens_ancestor_21.fa", "sites.tsv",
                      "anc_chr21.tsv.gz", chrom="chr21")
  rates = polarization_rates(counts)
  # CLI
  python3 epo_join.py --fasta ancestor_21.fa --sites sites.tsv \\
          --out anc_chr21.tsv.gz --chrom chr21 --summary counts.json
  gzip -dc calls.vcf.gz | grep -v '^##' | cut -f1,2,4,5 | \\
      python3 epo_join.py --fasta ancestor_21.fa --sites - --out anc.tsv.gz \\
          --chrom 21
Standard library only.
"""
import argparse
import gzip
import json
import os
import sys

HIGH = frozenset("ACGT")
LOW = frozenset("acgt")

OUT_HEADER = ("pos", "ref", "alt", "anc", "conf")

CONVENTION = ("EPO ancestral FASTA: uppercase ACGT = high-confidence call, "
              "lowercase acgt = low-confidence call, 'N'/'-'/'.' = no call")

__all__ = ["HIGH", "LOW", "OUT_HEADER", "CONVENTION", "load_fasta",
           "ancestral_base", "classify", "orient", "iter_sites",
           "join_sites", "polarization_rates", "main"]


# ------------------------------------------------------------------ FASTA
def load_fasta(path):
    """Read the FIRST record of a FASTA file.

    Returns (headers, seq): headers is the list of '>' lines seen (at most
    two: the first record's and, if present, the next one, at which reading
    stops), seq the first record's sequence with line breaks removed and
    case preserved. One chromosome per file is the intended use."""
    headers, chunks = [], []
    with open(path, "rt") as fh:
        cur = None
        for line in fh:
            if line.startswith(">"):
                headers.append(line.rstrip("\n"))
                if len(headers) > 1:
                    break  # only the first sequence is the chromosome
                cur = chunks
            elif cur is not None and len(headers) == 1:
                cur.append(line.strip())
    return headers, "".join(chunks)


def ancestral_base(seq, pos):
    """Ancestral base at 1-based position pos, or '.' when pos is outside
    [1, len(seq)]."""
    if 1 <= pos <= len(seq):
        return seq[pos - 1]
    return "."


# ------------------------------------------------------------- convention
def classify(anc):
    """EPO ancestral base -> 'high' | 'low' | 'none'."""
    if anc in HIGH:
        return "high"
    if anc in LOW:
        return "low"
    return "none"


def orient(ref, alt, anc):
    """Orient a site against its ancestral base.

    Returns (key, ancestral, derived):
      key = None            not a biallelic SNV (len != 1 or not in ACGT):
                            the site is written to the table but not
                            entered in the orientation counters;
      key = 'unpolarizable' anc is no call ('none' class);
      key = 'anc_eq_ref'    REF ancestral, ALT derived;
      key = 'anc_eq_alt'    ALT ancestral, REF derived;
      key = 'anc_neither'   the (case-folded) ancestral base is a third
                            allele: mismatch class, not polarizable.
    ancestral/derived are the oriented alleles for the two anc_eq_* keys and
    None otherwise. Case of anc is ignored for the comparison (a low-
    confidence 't' orients like 'T'); the confidence class is classify(anc).
    """
    if not (len(ref) == 1 and len(alt) == 1 and ref in HIGH and alt in HIGH):
        return None, None, None
    if classify(anc) == "none":
        return "unpolarizable", None, None
    au = anc.upper()
    if au == ref:
        return "anc_eq_ref", ref, alt
    if au == alt:
        return "anc_eq_alt", alt, ref
    return "anc_neither", None, None


# ------------------------------------------------------------------ sites
def _opener(path):
    if path == "-":
        return sys.stdin
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path, "rt")


def _is_int(s):
    try:
        int(s)
        return True
    except ValueError:
        return False


def _chrom_forms(chrom):
    """The accepted spellings of a chromosome name: with and without 'chr'."""
    c = str(chrom)
    bare = c[3:] if c.lower().startswith("chr") else c
    return {c, bare, "chr" + bare}


def iter_sites(source, layout="auto", chrom=None, counters=None):
    """Yield (pos_str, ref, alt) for every retained row of a sites table.

    source   path ('.gz' accepted, '-' = stdin) or an open text handle.
    layout   'auto' | 'chrom' | 'nochrom' (see module docstring).
    chrom    with layout 'chrom': keep only rows of this chromosome ('chrN'
             and 'N' both accepted); None keeps every row.
    counters optional dict; 'skipped_other_chrom' and 'skipped_short_rows'
             are incremented in place, and 'layout' records the layout used.
    Lines starting with '#' are skipped; a header row is consumed ('nochrom':
    first three names exactly pos ref alt; 'chrom': second name 'pos' in any
    case); rows with too few columns are skipped and counted."""
    if counters is None:
        counters = {}
    counters.setdefault("skipped_other_chrom", 0)
    counters.setdefault("skipped_short_rows", 0)
    keep = _chrom_forms(chrom) if chrom is not None else None
    fh = source if hasattr(source, "readline") else _opener(source)
    close = fh is not source and fh is not sys.stdin
    try:
        lay = layout
        first = True
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if first:
                first = False
                if lay == "auto":
                    if parts[0] == "pos":
                        lay = "nochrom"
                    elif (_is_int(parts[0]) and
                          (len(parts) < 4 or not _is_int(parts[1]))):
                        lay = "nochrom"
                    else:
                        lay = "chrom"
                counters["layout"] = lay
                if lay not in ("chrom", "nochrom"):
                    raise ValueError(f"unknown layout {lay!r}")
                if lay == "nochrom" and parts[0] == "pos":
                    if parts[:3] != ["pos", "ref", "alt"]:
                        raise ValueError(
                            f"sites header must start 'pos ref alt', got "
                            f"{parts[:3]!r}")
                    continue                      # header row consumed
                if (lay == "chrom" and len(parts) >= 4
                        and parts[1].lower() == "pos"):
                    continue                      # 'chrom pos ref alt' header
            if lay == "chrom":
                if len(parts) < 4:
                    counters["skipped_short_rows"] += 1
                    continue
                c, pos_s, ref, alt = parts[0], parts[1], parts[2], parts[3]
                if keep is not None and c not in keep:
                    counters["skipped_other_chrom"] += 1
                    continue
            else:
                if len(parts) < 3:
                    counters["skipped_short_rows"] += 1
                    continue
                pos_s, ref, alt = parts[0], parts[1], parts[2]
            yield pos_s, ref, alt
        if first:
            counters.setdefault("layout", layout if layout != "auto" else "chrom")
    finally:
        if close:
            fh.close()


# ------------------------------------------------------------------- join
def _new_register():
    return {"n_biallelic_snv": 0,
            "high": {"anc_eq_ref": 0, "anc_eq_alt": 0, "anc_neither": 0},
            "low": {"anc_eq_ref": 0, "anc_eq_alt": 0, "anc_neither": 0},
            "unpolarizable": 0}


def join_sites(fasta, sites, out=None, chrom=None, layout="auto",
               compresslevel=6, flush_rows=100000):
    """Join a sites table against an EPO ancestral FASTA and write the
    polarized table.

    fasta   per-chromosome ancestral-allele FASTA (first record used).
    sites   sites table: path ('.gz' ok, '-' = stdin), open text handle, or
            an iterable of (pos, ref, alt) / (chrom, pos, ref, alt) tuples.
    out     output path for the 'pos ref alt anc conf' table ('.gz' ->
            gzip at compresslevel, '-' -> stdout, None -> no table, counters
            only) or an open text handle.
    chrom   chromosome filter for 'chrom'-layout tables ('chrN' and 'N' both
            match); None keeps every row. Recorded in the counters.
    Returns the counters dict:
      chrom, layout, fasta_header, fasta_extra_sequences (1 if the FASTA
      holds records beyond the first, else 0), fasta_length, total_sites,
      out_of_fasta_range, skipped_other_chrom, skipped_short_rows,
      conf_counts_all_sites {high, low, none}, biallelic_snv {n_biallelic_snv,
      high{anc_eq_ref, anc_eq_alt, anc_neither}, low{...}, unpolarizable},
      mismatch_class_biallelic_snv_anc_neither, n_anc_called, convention.
    """
    headers, seq = load_fasta(fasta)
    L = len(seq)

    aux = {}
    if isinstance(sites, os.PathLike):
        sites = os.fspath(sites)
    if isinstance(sites, str) or hasattr(sites, "readline"):
        rows = iter_sites(sites, layout=layout, chrom=chrom, counters=aux)
    else:
        keep = _chrom_forms(chrom) if chrom is not None else None
        aux.update(skipped_other_chrom=0, skipped_short_rows=0,
                   layout="tuples")

        def rows_from_tuples(it):
            for v in it:
                v = tuple(v)
                if len(v) >= 4:
                    if keep is not None and str(v[0]) not in keep:
                        aux["skipped_other_chrom"] += 1
                        continue
                    yield str(v[1]), v[2], v[3]
                elif len(v) == 3:
                    yield str(v[0]), v[1], v[2]
                else:
                    aux["skipped_short_rows"] += 1
        rows = rows_from_tuples(sites)

    n_total = 0
    n_out_of_range = 0
    conf_counts = {"high": 0, "low": 0, "none": 0}
    b = _new_register()

    close_out = False
    if out is None:
        ofh = None
    elif hasattr(out, "write"):
        ofh = out
    elif out == "-":
        ofh = sys.stdout
    elif str(out).endswith(".gz"):
        ofh = gzip.open(out, "wt", compresslevel=compresslevel)
        close_out = True
    else:
        ofh = open(out, "wt")
        close_out = True
    try:
        if ofh is not None:
            ofh.write("\t".join(OUT_HEADER) + "\n")
        wbuf = []
        for pos_s, ref, alt in rows:
            n_total += 1
            pos = int(pos_s)
            if 1 <= pos <= L:
                anc = seq[pos - 1]
            else:
                anc = "."
                n_out_of_range += 1
            conf = classify(anc)
            conf_counts[conf] += 1
            if ofh is not None:
                wbuf.append(f"{pos}\t{ref}\t{alt}\t{anc}\t{conf}\n")
                if len(wbuf) >= flush_rows:
                    ofh.write("".join(wbuf))
                    wbuf = []
            key, _a, _d = orient(ref, alt, anc)
            if key is not None:
                b["n_biallelic_snv"] += 1
                if key == "unpolarizable":
                    b["unpolarizable"] += 1
                else:
                    b[conf][key] += 1
        if ofh is not None:
            ofh.write("".join(wbuf))
    finally:
        if close_out:
            ofh.close()

    n_neither = b["high"]["anc_neither"] + b["low"]["anc_neither"]
    return {
        "chrom": chrom,
        "layout": aux.get("layout"),
        "fasta_header": headers[0] if headers else None,
        "fasta_extra_sequences": max(0, len(headers) - 1),
        "fasta_length": L,
        "total_sites": n_total,
        "out_of_fasta_range": n_out_of_range,
        "skipped_other_chrom": aux.get("skipped_other_chrom", 0),
        "skipped_short_rows": aux.get("skipped_short_rows", 0),
        "conf_counts_all_sites": conf_counts,
        "biallelic_snv": b,
        "mismatch_class_biallelic_snv_anc_neither": n_neither,
        "n_anc_called": b["n_biallelic_snv"] - b["unpolarizable"],
        "convention": CONVENTION,
    }


def polarization_rates(counts):
    """Summary fractions of the biallelic-SNV orientation counters.

    Returns n_biallelic_snv, n_anc_called (= biallelic SNVs with a high or
    low call), frac_anc_eq_ref_of_called, frac_anc_eq_alt_of_called,
    frac_anc_neither_of_called (the ancestral-mismatch rate),
    frac_low_conf_of_called, frac_unpolarizable_of_biallelic, and the
    by_confidence sub-registers. A fraction with a zero denominator is
    None."""
    b = counts["biallelic_snv"]
    n_called = b["n_biallelic_snv"] - b["unpolarizable"]
    n_low = (b["low"]["anc_eq_ref"] + b["low"]["anc_eq_alt"]
             + b["low"]["anc_neither"])
    n_neither = b["high"]["anc_neither"] + b["low"]["anc_neither"]

    def frac(x, d):
        return (x / d) if d else None

    return {
        "n_biallelic_snv": b["n_biallelic_snv"],
        "n_anc_called": n_called,
        "frac_anc_eq_ref_of_called": frac(
            b["high"]["anc_eq_ref"] + b["low"]["anc_eq_ref"], n_called),
        "frac_anc_eq_alt_of_called": frac(
            b["high"]["anc_eq_alt"] + b["low"]["anc_eq_alt"], n_called),
        "frac_anc_neither_of_called": frac(n_neither, n_called),
        "frac_low_conf_of_called": frac(n_low, n_called),
        "frac_unpolarizable_of_biallelic": frac(
            b["unpolarizable"], b["n_biallelic_snv"]),
        "by_confidence": {"high": dict(b["high"]), "low": dict(b["low"])},
    }


# -------------------------------------------------------------------- CLI
def main(argv=None):
    ap = argparse.ArgumentParser(
        description="Join a sites table against an EPO ancestral-allele "
                    "FASTA; write 'pos ref alt anc conf' and print counters.")
    ap.add_argument("--fasta", required=True,
                    help="per-chromosome ancestral FASTA (first record used)")
    ap.add_argument("--sites", required=True,
                    help="sites TSV: chrom pos ref alt [...] or pos ref alt "
                         "[...] ('#' lines skipped; .gz ok; '-' = stdin)")
    ap.add_argument("--out", required=True,
                    help="output table (.gz -> gzip; '-' = stdout)")
    ap.add_argument("--chrom", default=None,
                    help="keep only this chromosome's rows in a chrom-layout "
                         "table ('chrN' and 'N' both match)")
    ap.add_argument("--layout", default="auto",
                    choices=("auto", "chrom", "nochrom"))
    ap.add_argument("--summary", default=None,
                    help="write the counters (+ rates) as JSON to this path")
    a = ap.parse_args(argv)

    counts = join_sites(a.fasta, a.sites, a.out, chrom=a.chrom,
                        layout=a.layout)
    rates = polarization_rates(counts)
    if a.summary:
        doc = dict(counts)
        doc["rates"] = rates
        doc["fasta"] = a.fasta
        doc["sites"] = a.sites
        doc["out"] = a.out
        with open(a.summary, "w") as fh:
            json.dump(doc, fh, indent=1)
            fh.write("\n")
    cc, b = counts["conf_counts_all_sites"], counts["biallelic_snv"]
    dest = sys.stderr if a.out == "-" else sys.stdout
    print(f"epo_join {counts['chrom'] or '*'}: {counts['total_sites']} sites; "
          f"high {cc['high']} low {cc['low']} none {cc['none']}; "
          f"out_of_range {counts['out_of_fasta_range']} skipped_other_chrom "
          f"{counts['skipped_other_chrom']}; biallelic {b['n_biallelic_snv']} "
          f"called {counts['n_anc_called']} neither "
          f"{counts['mismatch_class_biallelic_snv_anc_neither']} unpol "
          f"{b['unpolarizable']}", file=dest, flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
