/* posq_kernel.c -- multi-topology certified positive Gauss-Jacobi tensor
 * sweep kernel for the posq gate-quartet production run.
 *
 * REGISTER-PRESERVING: same mathematical chain as the proven python engine
 * (bench_rung1_gauss.sweep_range / rung-4b msweep): per-node VALUE evaluation
 * P_y = A_y + B_y*v (never coefficient expansion), positive-weighted ball
 * sum in arb at the job's precision. Only the arithmetic substrate changes
 * (python-flint -> C flint/arb) plus a shared squaring table for the
 * per-topology powers loop (exact same ball semantics, different op order --
 * validated by ball agreement + exact-rational surrogate + recorded-node
 * containment controls in the driver).
 *
 * Modes: cat  -- caterpillar ab_fiber (16 tree-slot pattern values)
 *        bal  -- balanced min/diff wedge fiber (S1,S2,T collapse)
 *        ms   -- majorant sweep, Bernstein-table fiber (topology-independent
 *                Ptilde values; works for cat and bal tables alike)
 * Multi-vector: NVEC exponent 16-vectors over the shared value table; one
 * certified ball per vector is emitted (man/exp pairs, exactly re-readable).
 *
 * Build: gcc -O2 -o posq_kernel posq_kernel.c -lflint -lmpfr -lgmp -lm
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <flint/flint.h>
#include <flint/fmpz.h>
#include <flint/arf.h>
#include <flint/mag.h>
#include <flint/arb.h>

#define MAXN 1024
#define NV 16

static void arb_from_manexp(arb_t x, const fmpz_t mm, const fmpz_t me,
                            const fmpz_t rm, const fmpz_t re)
{
    arf_t t;
    arf_init(t);
    arf_set_fmpz_2exp(arb_midref(x), mm, me);
    arf_set_fmpz_2exp(t, rm, re);
    arf_get_mag(arb_radref(x), t);      /* rounds UP: sound */
    arf_clear(t);
}

static void read_ball_line(FILE *f, arb_t x)
{
    fmpz_t mm, me, rm, re;
    fmpz_init(mm); fmpz_init(me); fmpz_init(rm); fmpz_init(re);
    if (fmpz_fread(f, mm) <= 0 || fmpz_fread(f, me) <= 0 ||
        fmpz_fread(f, rm) <= 0 || fmpz_fread(f, re) <= 0)
    { fprintf(stderr, "KERNEL: ball read failed\n"); exit(2); }
    arb_from_manexp(x, mm, me, rm, re);
    fmpz_clear(mm); fmpz_clear(me); fmpz_clear(rm); fmpz_clear(re);
}

static void write_ball(FILE *f, const char *tag, int idx, const arb_t x)
{
    fmpz_t mm, me, rm, re;
    arf_t t;
    fmpz_init(mm); fmpz_init(me); fmpz_init(rm); fmpz_init(re);
    arf_init(t);
    if (arf_is_zero(arb_midref(x))) { fmpz_zero(mm); fmpz_zero(me); }
    else arf_get_fmpz_2exp(mm, me, arb_midref(x));
    arf_set_mag(t, arb_radref(x));
    if (arf_is_zero(t)) { fmpz_zero(rm); fmpz_zero(re); }
    else arf_get_fmpz_2exp(rm, re, t);
    char *s1 = fmpz_get_str(NULL, 10, mm);
    char *s2 = fmpz_get_str(NULL, 10, me);
    char *s3 = fmpz_get_str(NULL, 10, rm);
    char *s4 = fmpz_get_str(NULL, 10, re);
    fprintf(f, "%s %d %s %s %s %s\n", tag, idx, s1, s2, s3, s4);
    flint_free(s1); flint_free(s2); flint_free(s3); flint_free(s4);
    arf_clear(t);
    fmpz_clear(mm); fmpz_clear(me); fmpz_clear(rm); fmpz_clear(re);
}

typedef struct { int n; arb_ptr x; arb_ptr w; } rule_t;

static void read_rule(FILE *f, rule_t *r)
{
    if (fscanf(f, "%d", &r->n) != 1 || r->n < 1 || r->n > MAXN)
    { fprintf(stderr, "KERNEL: bad rule n\n"); exit(2); }
    r->x = _arb_vec_init(r->n);
    r->w = _arb_vec_init(r->n);
    for (int i = 0; i < r->n; i++) read_ball_line(f, r->x + i);
    for (int i = 0; i < r->n; i++) read_ball_line(f, r->w + i);
}

int main(int argc, char **argv)
{
    if (argc != 2) { fprintf(stderr, "usage: posq_kernel <jobfile>\n"); return 2; }
    FILE *jf = fopen(argv[1], "r");
    if (!jf) { fprintf(stderr, "KERNEL: no jobfile\n"); return 2; }
    char key[64], mode[16], rulepath[4096], tabpath[4096], outpath[4096];
    long prec = 0; long i0 = -1, i1 = -1; int nvec = 0;
    fmpz_t pnum, pden; fmpz_init(pnum); fmpz_init(pden);
    int vec[64][NV]; memset(vec, 0, sizeof vec);
    mode[0] = tabpath[0] = rulepath[0] = outpath[0] = 0;
    while (fscanf(jf, "%63s", key) == 1) {
        if (!strcmp(key, "PREC")) { if (fscanf(jf, "%ld", &prec) != 1) return 2; }
        else if (!strcmp(key, "MODE")) { if (fscanf(jf, "%15s", mode) != 1) return 2; }
        else if (!strcmp(key, "PNUM")) { if (fmpz_fread(jf, pnum) <= 0) return 2; }
        else if (!strcmp(key, "PDEN")) { if (fmpz_fread(jf, pden) <= 0) return 2; }
        else if (!strcmp(key, "I0")) { if (fscanf(jf, "%ld", &i0) != 1) return 2; }
        else if (!strcmp(key, "I1")) { if (fscanf(jf, "%ld", &i1) != 1) return 2; }
        else if (!strcmp(key, "RULES")) { if (fscanf(jf, "%4095s", rulepath) != 1) return 2; }
        else if (!strcmp(key, "TABLES")) { if (fscanf(jf, "%4095s", tabpath) != 1) return 2; }
        else if (!strcmp(key, "OUT")) { if (fscanf(jf, "%4095s", outpath) != 1) return 2; }
        else if (!strcmp(key, "NVEC")) {
            if (fscanf(jf, "%d", &nvec) != 1 || nvec < 1 || nvec > 64) return 2;
            for (int s = 0; s < nvec; s++)
                for (int m = 0; m < NV; m++)
                    if (fscanf(jf, "%d", &vec[s][m]) != 1) return 2;
        }
        else { fprintf(stderr, "KERNEL: bad key %s\n", key); return 2; }
    }
    fclose(jf);
    if (!prec || !mode[0] || !rulepath[0] || !outpath[0] || nvec < 1 || i0 < 0)
    { fprintf(stderr, "KERNEL: incomplete job\n"); return 2; }
    int is_cat = !strcmp(mode, "cat"), is_bal = !strcmp(mode, "bal"),
        is_ms = !strcmp(mode, "ms");
    if (!is_cat && !is_bal && !is_ms)
    { fprintf(stderr, "KERNEL: bad mode\n"); return 2; }

    FILE *rf = fopen(rulepath, "r");
    if (!rf) { fprintf(stderr, "KERNEL: no rules file\n"); return 2; }
    rule_t R1, R2, R3;
    read_rule(rf, &R1); read_rule(rf, &R2); read_rule(rf, &R3);
    fclose(rf);
    if (i1 > R1.n) { fprintf(stderr, "KERNEL: i1 > n1\n"); return 2; }

    /* pi1 = pnum/pden ball at prec (contains exact rational; sweep semantics
       identical to python arb(P0_NUM)/P0_DEN) -- unused in ms mode */
    arb_t pi0, pi1;
    arb_init(pi0); arb_init(pi1);
    arb_set_fmpz(pi1, pnum);
    arb_div_fmpz(pi1, pi1, pden, prec);
    arb_one(pi0); arb_sub(pi0, pi0, pi1, prec);

    /* Bernstein tables for ms mode: 16 values x 5 x 4 x {l=0,1} rationals */
    static arb_t D0[NV][5][4], D1[NV][5][4];
    if (is_ms) {
        FILE *tf = fopen(tabpath, "r");
        if (!tf) { fprintf(stderr, "KERNEL: no tables file\n"); return 2; }
        fmpz_t nu, de; fmpz_init(nu); fmpz_init(de);
        for (int y = 0; y < NV; y++)
            for (int i = 0; i < 5; i++)
                for (int j = 0; j < 4; j++)
                    for (int l = 0; l < 2; l++) {
                        if (fmpz_fread(tf, nu) <= 0 || fmpz_fread(tf, de) <= 0)
                        { fprintf(stderr, "KERNEL: table read fail\n"); return 2; }
                        arb_ptr d = l ? D1[y][i][j] : D0[y][i][j];
                        arb_init(d);
                        arb_set_fmpz(d, nu);
                        arb_div_fmpz(d, d, de, prec);
                    }
        fmpz_clear(nu); fmpz_clear(de);
        fclose(tf);
    }

    /* used value indices + squaring levels + per-vector plans */
    int used[NV], nused = 0, uidx[NV];
    for (int m = 0; m < NV; m++) uidx[m] = -1;
    for (int m = 0; m < NV; m++) {
        int hit = 0;
        for (int s = 0; s < nvec; s++) if (vec[s][m]) hit = 1;
        if (hit) { uidx[m] = nused; used[nused++] = m; }
    }
    int maxlev[NV];
    for (int u = 0; u < nused; u++) {
        int lv = 0;
        for (int s = 0; s < nvec; s++) {
            int n = vec[s][used[u]];
            int b = 0; while ((1 << (b + 1)) <= n) b++;
            if (n && b > lv) lv = b;
        }
        maxlev[u] = lv;
    }
    /* plan[s]: flattened (u, lev) pairs over set bits of each exponent */
    int plan[64][NV * 8][2], plen[64];
    for (int s = 0; s < nvec; s++) {
        plen[s] = 0;
        for (int u = 0; u < nused; u++) {
            int n = vec[s][used[u]];
            for (int b = 0; b <= maxlev[u]; b++)
                if (n & (1 << b)) {
                    plan[s][plen[s]][0] = u;
                    plan[s][plen[s]][1] = b;
                    plen[s]++;
                }
        }
    }

    /* working vars */
    arb_t one, t1, t2, e1, e2, m1, m2, u12, m12, F;
    arb_init(one); arb_one(one);
    arb_init(t1); arb_init(t2); arb_init(e1); arb_init(e2);
    arb_init(m1); arb_init(m2); arb_init(u12); arb_init(m12); arb_init(F);
    arb_t P1[2][2], P2[2][2], P12[2][2], gg[2][2][2];
    for (int a = 0; a < 2; a++) for (int b = 0; b < 2; b++) {
        arb_init(P1[a][b]); arb_init(P2[a][b]); arb_init(P12[a][b]);
        for (int c = 0; c < 2; c++) arb_init(gg[a][b][c]);
    }
    arb_t C1a[4][2], C2a[4][2], S1[4], S2[4];
    for (int q = 0; q < 4; q++) {
        arb_init(S1[q]); arb_init(S2[q]);
        for (int a = 0; a < 2; a++) { arb_init(C1a[q][a]); arb_init(C2a[q][a]); }
    }
    arb_t B1b[5], B2b[4];
    for (int i = 0; i < 5; i++) arb_init(B1b[i]);
    for (int j = 0; j < 4; j++) arb_init(B2b[j]);
    static const int CB1[5] = {1, 4, 6, 4, 1};
    static const int CB2[4] = {1, 3, 3, 1};
    arb_ptr A = _arb_vec_init(nused), Bv = _arb_vec_init(nused);
    arb_t S[NV][8];
    for (int u = 0; u < nused; u++)
        for (int b = 0; b <= maxlev[u]; b++) arb_init(S[u][b]);
    arb_ptr tot = _arb_vec_init(nvec), row = _arb_vec_init(nvec),
            fib = _arb_vec_init(nvec);

    clock_t ck0 = clock();
    long nodes_done = 0;
    for (long i = i0; i < i1; i++) {
        arb_srcptr u1 = R1.x + i;
        _arb_vec_zero(row, nvec);
        for (long j = 0; j < R2.n; j++) {
            arb_srcptr u2 = R2.x + j;
            /* ---------------- fiber stage ---------------- */
            if (is_cat) {
                arb_sub(m1, one, u1, prec);
                arb_sub(m2, one, u2, prec);
                arb_mul(u12, u1, u2, prec);
                arb_sub(m12, one, u12, prec);
                /* P[i][j] rows for decay e: ((pi0+pi1*e, pi1*(1-e)),
                                             (pi0*(1-e), pi1+pi0*e)) */
                arb_mul(P1[0][0], pi1, u1, prec); arb_add(P1[0][0], P1[0][0], pi0, prec);
                arb_mul(P1[0][1], pi1, m1, prec);
                arb_mul(P1[1][0], pi0, m1, prec);
                arb_mul(P1[1][1], pi0, u1, prec); arb_add(P1[1][1], P1[1][1], pi1, prec);
                arb_mul(P2[0][0], pi1, u2, prec); arb_add(P2[0][0], P2[0][0], pi0, prec);
                arb_mul(P2[0][1], pi1, m2, prec);
                arb_mul(P2[1][0], pi0, m2, prec);
                arb_mul(P2[1][1], pi0, u2, prec); arb_add(P2[1][1], P2[1][1], pi1, prec);
                arb_mul(P12[0][0], pi1, u12, prec); arb_add(P12[0][0], P12[0][0], pi0, prec);
                arb_mul(P12[0][1], pi1, m12, prec);
                arb_mul(P12[1][0], pi0, m12, prec);
                arb_mul(P12[1][1], pi0, u12, prec); arb_add(P12[1][1], P12[1][1], pi1, prec);
                for (int w = 0; w < 2; w++)
                    for (int a = 0; a < 2; a++)
                        for (int b = 0; b < 2; b++) {
                            arb_mul(t1, P1[0][a], P1[0][b], prec);
                            arb_mul(t1, t1, P2[w][0], prec);
                            arb_mul(t2, P1[1][a], P1[1][b], prec);
                            arb_mul(t2, t2, P2[w][1], prec);
                            arb_add(gg[w][a][b], t1, t2, prec);
                        }
                for (int u = 0; u < nused; u++) {
                    int y = used[u];
                    int a = (y >> 3) & 1, b = (y >> 2) & 1,
                        c = (y >> 1) & 1, d = y & 1;
                    /* base0/base1 */
                    arb_mul(t1, gg[0][a][b], P12[0][c], prec);
                    arb_mul(t1, t1, pi0, prec);
                    arb_mul(t2, gg[1][a][b], P12[1][c], prec);
                    arb_mul(t2, t2, pi1, prec);
                    arb_add(A + u, t1, t2, prec);
                    if (d == 0) {
                        /* Bv = (b0*pi1 - b1*pi0) * u12 ; A *= pi0 */
                        arb_mul(Bv + u, t1, pi1, prec);
                        arb_mul(F, t2, pi0, prec);
                        arb_sub(Bv + u, Bv + u, F, prec);
                        arb_mul(A + u, A + u, pi0, prec);
                    } else {
                        /* Bv = (b1*pi0 - b0*pi1) * u12 ; A *= pi1 */
                        arb_mul(Bv + u, t2, pi0, prec);
                        arb_mul(F, t1, pi1, prec);
                        arb_sub(Bv + u, Bv + u, F, prec);
                        arb_mul(A + u, A + u, pi1, prec);
                    }
                    arb_mul(Bv + u, Bv + u, u12, prec);
                }
            } else if (is_bal) {
                /* e1 = w1*w2 (old leaves), e2 = w1 (young leaves) */
                arb_mul(e1, u1, u2, prec);
                arb_set(e2, u1);
                arb_sub(m1, one, e1, prec);
                arb_sub(m2, one, e2, prec);
                arb_mul(P1[0][0], pi1, e1, prec); arb_add(P1[0][0], P1[0][0], pi0, prec);
                arb_mul(P1[0][1], pi1, m1, prec);
                arb_mul(P1[1][0], pi0, m1, prec);
                arb_mul(P1[1][1], pi0, e1, prec); arb_add(P1[1][1], P1[1][1], pi1, prec);
                arb_mul(P2[0][0], pi1, e2, prec); arb_add(P2[0][0], P2[0][0], pi0, prec);
                arb_mul(P2[0][1], pi1, m2, prec);
                arb_mul(P2[1][0], pi0, m2, prec);
                arb_mul(P2[1][1], pi0, e2, prec); arb_add(P2[1][1], P2[1][1], pi1, prec);
                for (int xa = 0; xa < 2; xa++)
                    for (int xb = 0; xb < 2; xb++) {
                        int q = 2 * xa + xb;
                        arb_mul(C1a[q][0], P1[0][xa], P1[0][xb], prec);
                        arb_mul(C1a[q][1], P1[1][xa], P1[1][xb], prec);
                        arb_mul(C2a[q][0], P2[0][xa], P2[0][xb], prec);
                        arb_mul(C2a[q][1], P2[1][xa], P2[1][xb], prec);
                        arb_mul(t1, C1a[q][0], pi0, prec);
                        arb_mul(t2, C1a[q][1], pi1, prec);
                        arb_add(S1[q], t1, t2, prec);
                        arb_mul(t1, C2a[q][0], pi0, prec);
                        arb_mul(t2, C2a[q][1], pi1, prec);
                        arb_add(S2[q], t1, t2, prec);
                    }
                for (int u = 0; u < nused; u++) {
                    int y = used[u];
                    int q12 = (y >> 2) & 3, q34 = y & 3;
                    arb_mul(A + u, S1[q12], S2[q34], prec);
                    /* T = pi0*C1[0]*C2[0] + pi1*C1[1]*C2[1] */
                    arb_mul(t1, C1a[q12][0], C2a[q34][0], prec);
                    arb_mul(t1, t1, pi0, prec);
                    arb_mul(t2, C1a[q12][1], C2a[q34][1], prec);
                    arb_mul(t2, t2, pi1, prec);
                    arb_add(t1, t1, t2, prec);
                    arb_sub(t1, t1, A + u, prec);
                    arb_mul(Bv + u, t1, u2, prec);   /* * w2 */
                }
            } else { /* ms: Bernstein-table fiber */
                arb_sub(m1, one, u1, prec);
                arb_sub(m2, one, u2, prec);
                for (int a = 0; a <= 4; a++) {
                    arb_pow_ui(t1, u1, a, prec);
                    arb_pow_ui(t2, m1, 4 - a, prec);
                    arb_mul(B1b[a], t1, t2, prec);
                    arb_mul_ui(B1b[a], B1b[a], CB1[a], prec);
                }
                for (int b = 0; b <= 3; b++) {
                    arb_pow_ui(t1, u2, b, prec);
                    arb_pow_ui(t2, m2, 3 - b, prec);
                    arb_mul(B2b[b], t1, t2, prec);
                    arb_mul_ui(B2b[b], B2b[b], CB2[b], prec);
                }
                for (int u = 0; u < nused; u++) {
                    int y = used[u];
                    arb_zero(A + u); arb_zero(Bv + u);
                    for (int a = 0; a <= 4; a++)
                        for (int b = 0; b <= 3; b++) {
                            arb_mul(t1, B1b[a], B2b[b], prec);
                            arb_addmul(A + u, D0[y][a][b], t1, prec);
                            arb_addmul(Bv + u, D1[y][a][b], t1, prec);
                        }
                    arb_sub(Bv + u, Bv + u, A + u, prec);  /* Bt = A1 - A0 */
                }
            }
            /* ---------------- v-fiber + powers ---------------- */
            _arb_vec_zero(fib, nvec);
            for (long k = 0; k < R3.n; k++) {
                arb_srcptr v = R3.x + k;
                for (int u = 0; u < nused; u++) {
                    arb_mul(S[u][0], Bv + u, v, prec);
                    arb_add(S[u][0], S[u][0], A + u, prec);
                    for (int b = 1; b <= maxlev[u]; b++)
                        arb_sqr(S[u][b], S[u][b - 1], prec);
                }
                for (int s = 0; s < nvec; s++) {
                    arb_set(F, S[plan[s][0][0]][plan[s][0][1]]);
                    for (int q = 1; q < plen[s]; q++)
                        arb_mul(F, F, S[plan[s][q][0]][plan[s][q][1]], prec);
                    arb_addmul(fib + s, R3.w + k, F, prec);
                }
                nodes_done++;
            }
            for (int s = 0; s < nvec; s++)
                arb_addmul(row + s, R2.w + j, fib + s, prec);
        }
        for (int s = 0; s < nvec; s++)
            arb_addmul(tot + s, R1.w + i, row + s, prec);
    }
    double cpu_s = (double)(clock() - ck0) / CLOCKS_PER_SEC;

    FILE *of = fopen(outpath, "w");
    if (!of) { fprintf(stderr, "KERNEL: cannot open out\n"); return 2; }
    fprintf(of, "MODE %s PREC %ld I0 %ld I1 %ld NVEC %d NODES %ld CPU_S %.1f\n",
            mode, prec, i0, i1, nvec, nodes_done, cpu_s);
    for (int s = 0; s < nvec; s++)
        write_ball(of, "TOT", s, tot + s);
    fclose(of);
    printf("[posq_kernel %s] nvec=%d nodes=%ld cpu=%.1fs %.3fus/node\n",
           mode, nvec, nodes_done, cpu_s, 1e6 * cpu_s / (nodes_done ? nodes_done : 1));
    return 0;
}
