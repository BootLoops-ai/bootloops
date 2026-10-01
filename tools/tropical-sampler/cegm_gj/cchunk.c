/* cchunk.c — C-MPFR inner loop for the X(3,6) cone evaluator.
 *
 * Implements EXACTLY the hoisted algorithm of t3_prod.py:eval_chunk_hoist:
 * same Gauss-Jacobi nodes/weights (computed by the certified Python code and
 * passed in as decimal strings), same factor tables, same summation
 * structure.  Support-hoisting: factors depending on fewer than all 4 axes
 * get gam*log(Q) precomputed on their own subgrid; full-support factors pay
 * a live log per point.  Gate before production: must reproduce the Python
 * path at the roundoff floor.
 *
 * build: gcc -O2 -shared -fPIC -o cchunk.so cchunk.c -lmpfr -lgmp
 */
#include <mpfr.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

/* value = product over axes of powt[ax][idx[ax]][e[ax]] summed over monos,
 * then gam*log(Q).  powt built at prec. */

typedef struct {
    mpfr_t gam;
    int nm;
    int *e;        /* nm x 4 */
    int mask;      /* OR of axes with any e>0 */
} Fac;

static void qlog(mpfr_t out, Fac *f, mpfr_t **powt[4],
                 int i0, int i1, int i2, int i3, int prec)
{
    int idx[4] = {i0, i1, i2, i3};
    mpfr_t Q, t;
    mpfr_init2(Q, prec); mpfr_init2(t, prec);
    mpfr_set_ui(Q, 0, MPFR_RNDN);
    for (int m = 0; m < f->nm; m++) {
        mpfr_set_ui(t, 1, MPFR_RNDN);
        for (int ax = 0; ax < 4; ax++) {
            int e = f->e[4*m + ax];
            if (e)
                mpfr_mul(t, t, powt[ax][idx[ax]][e], MPFR_RNDN);
        }
        mpfr_add(Q, Q, t, MPFR_RNDN);
    }
    mpfr_log(Q, Q, MPFR_RNDN);
    mpfr_mul(out, Q, f->gam, MPFR_RNDN);
    mpfr_clear(Q); mpfr_clear(t);
}

/* Arguments (all strings decimal, exact round-trip digits):
 *   n, prec, lo, hi, det
 *   nodes: 4*n strings; wts: 4*n strings
 *   nfac, gam[nfac] strings "num/den", nm[nfac], monos flattened ints
 * Returns malloc'd decimal string of the chunk sum (caller frees). */
char *eval_chunk_c(int n, int prec, int lo, int hi, long det,
                   const char **nodes_s, const char **wts_s,
                   int nfac, const char **gam_s,
                   const int *nm_arr, const int *monos)
{
    mpfr_t nodes[4][n], wts[4][n];
    for (int ax = 0; ax < 4; ax++)
        for (int i = 0; i < n; i++) {
            mpfr_init2(nodes[ax][i], prec);
            mpfr_set_str(nodes[ax][i], nodes_s[ax*n + i], 10, MPFR_RNDN);
            mpfr_init2(wts[ax][i], prec);
            mpfr_set_str(wts[ax][i], wts_s[ax*n + i], 10, MPFR_RNDN);
        }
    /* factors */
    Fac *fac = malloc(nfac * sizeof(Fac));
    int off = 0;
    int maxe[4] = {0, 0, 0, 0};
    for (int f = 0; f < nfac; f++) {
        fac[f].nm = nm_arr[f];
        fac[f].e = (int *)(monos + off);
        off += 4 * nm_arr[f];
        fac[f].mask = 0;
        for (int m = 0; m < fac[f].nm; m++)
            for (int ax = 0; ax < 4; ax++)
                if (fac[f].e[4*m + ax]) {
                    fac[f].mask |= (1 << ax);
                    if (fac[f].e[4*m + ax] > maxe[ax])
                        maxe[ax] = fac[f].e[4*m + ax];
                }
        mpfr_init2(fac[f].gam, prec);
        mpq_t q; mpq_init(q);
        mpq_set_str(q, gam_s[f], 10); mpq_canonicalize(q);
        mpfr_set_q(fac[f].gam, q, MPFR_RNDN);
        mpq_clear(q);
    }
    /* power tables powt[ax][i][p], p = 0..maxe[ax] */
    mpfr_t **powt[4];
    for (int ax = 0; ax < 4; ax++) {
        powt[ax] = malloc(n * sizeof(mpfr_t *));
        for (int i = 0; i < n; i++) {
            powt[ax][i] = malloc((maxe[ax] + 1) * sizeof(mpfr_t));
            for (int p = 0; p <= maxe[ax]; p++) {
                mpfr_init2(powt[ax][i][p], prec);
                if (p == 0)
                    mpfr_set_ui(powt[ax][i][p], 1, MPFR_RNDN);
                else
                    mpfr_mul(powt[ax][i][p], powt[ax][i][p-1],
                             nodes[ax][i], MPFR_RNDN);
            }
        }
    }
    int nch = hi - lo;
    /* subgrid tables: for each factor with mask != 15, tabulate gam*log Q
     * over its axes.  Axis-0 index range is [lo,hi) (chunk), others full n.
     * Table layout: flat array, strides over the axes present in mask,
     * order axis0,1,2,3; axis-3 (if present) fastest. */
    mpfr_t *tab[nfac];
    long tsz[nfac];
    for (int f = 0; f < nfac; f++) {
        tab[f] = NULL;
        if (fac[f].mask == 15) continue;
        long sz = 1;
        for (int ax = 0; ax < 4; ax++)
            if (fac[f].mask & (1 << ax))
                sz *= (ax == 0 ? nch : n);
        tsz[f] = sz;
        tab[f] = malloc(sz * sizeof(mpfr_t));
        long c = 0;
        int r0 = (fac[f].mask & 1) ? nch : 1;
        int r1 = (fac[f].mask & 2) ? n : 1;
        int r2 = (fac[f].mask & 4) ? n : 1;
        int r3 = (fac[f].mask & 8) ? n : 1;
        for (int i0 = 0; i0 < r0; i0++)
        for (int i1 = 0; i1 < r1; i1++)
        for (int i2 = 0; i2 < r2; i2++)
        for (int i3 = 0; i3 < r3; i3++) {
            mpfr_init2(tab[f][c], prec);
            qlog(tab[f][c], &fac[f], powt,
                 (fac[f].mask & 1) ? lo + i0 : 0,
                 (fac[f].mask & 2) ? i1 : 0,
                 (fac[f].mask & 4) ? i2 : 0,
                 (fac[f].mask & 8) ? i3 : 0, prec);
            c++;
        }
    }
    /* main loop */
    mpfr_t tot, E0, E1, E2, E, acc, Q, w01, w012, tmp;
    mpfr_inits2(prec, tot, E0, E1, E2, E, acc, Q, w01, w012, tmp,
                (mpfr_ptr) 0);
    mpfr_set_ui(tot, 0, MPFR_RNDN);
    /* helper: table offset for factor f at (i0,i1,i2,i3) */
    #define TOFF(f, I0, I1, I2, I3) ({                                   \
        long _o = 0;                                                     \
        if (fac[f].mask & 1) _o = _o * nch + ((I0) - lo);                \
        if (fac[f].mask & 2) _o = _o * n + (I1);                         \
        if (fac[f].mask & 4) _o = _o * n + (I2);                         \
        if (fac[f].mask & 8) _o = _o * n + (I3);                         \
        _o; })
    for (int i0 = lo; i0 < hi; i0++) {
        mpfr_set_ui(E0, 0, MPFR_RNDN);
        for (int f = 0; f < nfac; f++)
            if (fac[f].mask != 15 && !(fac[f].mask & ~1))
                mpfr_add(E0, E0, tab[f][TOFF(f, i0, 0, 0, 0)], MPFR_RNDN);
        for (int i1 = 0; i1 < n; i1++) {
            mpfr_mul(w01, wts[0][i0], wts[1][i1], MPFR_RNDN);
            mpfr_set(E1, E0, MPFR_RNDN);
            for (int f = 0; f < nfac; f++)
                if (fac[f].mask != 15 && (fac[f].mask & 2)
                    && !(fac[f].mask & ~3))
                    mpfr_add(E1, E1, tab[f][TOFF(f, i0, i1, 0, 0)],
                             MPFR_RNDN);
            for (int i2 = 0; i2 < n; i2++) {
                mpfr_mul(w012, w01, wts[2][i2], MPFR_RNDN);
                mpfr_set(E2, E1, MPFR_RNDN);
                for (int f = 0; f < nfac; f++)
                    if (fac[f].mask != 15 && (fac[f].mask & 4)
                        && !(fac[f].mask & ~7))
                        mpfr_add(E2, E2, tab[f][TOFF(f, i0, i1, i2, 0)],
                                 MPFR_RNDN);
                mpfr_set_ui(acc, 0, MPFR_RNDN);
                for (int i3 = 0; i3 < n; i3++) {
                    mpfr_set(E, E2, MPFR_RNDN);
                    for (int f = 0; f < nfac; f++) {
                        if (fac[f].mask == 15) {
                            qlog(tmp, &fac[f], powt, i0, i1, i2, i3, prec);
                            mpfr_add(E, E, tmp, MPFR_RNDN);
                        } else if (fac[f].mask & 8) {
                            mpfr_add(E, E,
                                     tab[f][TOFF(f, i0, i1, i2, i3)],
                                     MPFR_RNDN);
                        }
                    }
                    mpfr_exp(E, E, MPFR_RNDN);
                    mpfr_mul(E, E, wts[3][i3], MPFR_RNDN);
                    mpfr_add(acc, acc, E, MPFR_RNDN);
                }
                mpfr_mul(acc, acc, w012, MPFR_RNDN);
                mpfr_add(tot, tot, acc, MPFR_RNDN);
            }
        }
    }
    mpfr_mul_si(tot, tot, det, MPFR_RNDN);
    /* output: enough digits for exact-ish round trip */
    int nd = (int)(prec * 0.30103) + 10;
    char fmt[32];
    snprintf(fmt, sizeof fmt, "%%.%dRe", nd);
    char *out = malloc(nd + 32);
    mpfr_sprintf(out, fmt, tot);
    /* cleanup */
    for (int f = 0; f < nfac; f++) {
        if (tab[f]) {
            for (long c = 0; c < tsz[f]; c++) mpfr_clear(tab[f][c]);
            free(tab[f]);
        }
        mpfr_clear(fac[f].gam);
    }
    free(fac);
    for (int ax = 0; ax < 4; ax++) {
        for (int i = 0; i < n; i++) {
            for (int p = 0; p <= maxe[ax]; p++)
                mpfr_clear(powt[ax][i][p]);
            free(powt[ax][i]);
            mpfr_clear(nodes[ax][i]); mpfr_clear(wts[ax][i]);
        }
        free(powt[ax]);
    }
    mpfr_clears(tot, E0, E1, E2, E, acc, Q, w01, w012, tmp, (mpfr_ptr) 0);
    mpfr_free_cache();
    return out;
}

void free_str(char *p) { free(p); }
