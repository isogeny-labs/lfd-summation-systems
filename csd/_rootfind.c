/* Point decomposition on y^2 + xy = x^3 + a x^2 + b over F_(2^n), n <= 63, without a Groebner basis.
 *
 * Kernels for csd/ic_rootfind.py (which holds the mathematics and builds every table used here):
 *   rf_fb     factor base {P : x(P) in V}
 *   rf_lrf    "linearized root finding": enumerate m-2 factor-base points, then solve the remaining
 *             two-summand problem as one F_2-linear system and one quadratic
 *   rf_brute  baseline: enumerate m-1 factor-base points, test x of the last point for membership in V
 *   rf_mitm4  baseline for m = 4: sorted table of pair sums, lookups for the other pair
 *
 * Linear maps F_2^n -> F_2^w are byte tables (8 x 256 words) built from the images of the basis bits.
 */
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <immintrin.h>

typedef uint64_t u64;
typedef struct { u64 t[8][256]; } lmap;

typedef struct {
    int n, k, w;            /* field degree, dim V, width of the linear system (n - dim V2 + 1) */
    u64 red, mask;          /* modulus minus x^n, 2^n - 1 */
    u64 a, b, sqrtb;
    u64 trmask;             /* Tr(x) = parity(x & trmask) */
    u64 v[64];              /* basis of V */
    u64 wv[64];             /* W(v_i), W = Q + id */
    u64 trv;                /* bit i = Tr(v_i) */
    lmap Q;                 /* Q(c)^2 + Q(c) = c for Tr(c) = 0 */
    lmap sqrt;              /* square root */
    lmap chk;               /* x in V  <=>  chk(x) == 0 */
    lmap pi;                /* F_2^n -> F_2^n / V2 */
} rf_ctx;

typedef struct {
    u64 targets, lin_solvable, lin_candidates, quad_split, found, degenerate, max_kernel_dim;
} rf_stats;

static inline u64 lm(const lmap *m, u64 x)
{
    return m->t[0][x & 255] ^ m->t[1][(x >> 8) & 255] ^ m->t[2][(x >> 16) & 255] ^ m->t[3][(x >> 24) & 255]
         ^ m->t[4][(x >> 32) & 255] ^ m->t[5][(x >> 40) & 255] ^ m->t[6][(x >> 48) & 255] ^ m->t[7][x >> 56];
}

static void lm_build(lmap *m, const u64 *images, int n)
{
    for (int j = 0; j < 8; j++) {
        m->t[j][0] = 0;
        for (int x = 1; x < 256; x++) {
            int low = __builtin_ctz(x), bit = 8 * j + low;
            m->t[j][x] = m->t[j][x & (x - 1)] ^ (bit < n ? images[bit] : 0);
        }
    }
}

static inline u64 fmul(const rf_ctx *c, u64 x, u64 y)
{
    __m128i p = _mm_clmulepi64_si128(_mm_cvtsi64_si128((long long)x), _mm_cvtsi64_si128((long long)y), 0);
    u64 lo = (u64)_mm_cvtsi128_si64(p), hi = (u64)_mm_extract_epi64(p, 1);
    u64 h = (hi << (64 - c->n)) | (lo >> c->n);
    lo &= c->mask;
    while (h) {
        __m128i q = _mm_clmulepi64_si128(_mm_cvtsi64_si128((long long)h), _mm_cvtsi64_si128((long long)c->red), 0);
        u64 ql = (u64)_mm_cvtsi128_si64(q), qh = (u64)_mm_extract_epi64(q, 1);
        lo ^= ql & c->mask;
        h = (qh << (64 - c->n)) | (ql >> c->n);
    }
    return lo;
}

static inline u64 fsq(const rf_ctx *c, u64 x) { return fmul(c, x, x); }

/* x^(2^n - 2) by Itoh-Tsujii: t_j = x^(2^j - 1), t_(i+j) = t_i^(2^j) t_j */
static u64 finv(const rf_ctx *c, u64 x)
{
    int e = c->n - 1, top = 31 - __builtin_clz((unsigned)e), have = 1;
    u64 t = x;
    for (int i = top - 1; i >= 0; i--) {
        u64 s = t;
        for (int j = 0; j < have; j++) s = fsq(c, s);
        t = fmul(c, s, t);
        have *= 2;
        if ((e >> i) & 1) { t = fmul(c, fsq(c, t), x); have++; }
    }
    return fsq(c, t);
}

static inline int ftr(const rf_ctx *c, u64 x) { return __builtin_parityll(x & c->trmask); }

typedef struct { u64 x, y; int inf; } pt;

static pt padd(const rf_ctx *c, pt P, pt Q)
{
    pt R = {0, 0, 0};
    if (P.inf) return Q;
    if (Q.inf) return P;
    u64 lam;
    if (P.x == Q.x) {
        if (Q.y == (P.x ^ P.y)) { R.inf = 1; return R; }
        lam = P.x ^ fmul(c, P.y, finv(c, P.x));
        R.x = fsq(c, lam) ^ lam ^ c->a;
        R.y = fsq(c, P.x) ^ fmul(c, lam ^ 1, R.x);
        return R;
    }
    lam = fmul(c, P.y ^ Q.y, finv(c, P.x ^ Q.x));
    R.x = fsq(c, lam) ^ lam ^ P.x ^ Q.x ^ c->a;
    R.y = fmul(c, lam, P.x ^ R.x) ^ R.x ^ P.y;
    return R;
}

static inline pt pneg(pt P) { P.y ^= P.x; return P; }

void rf_init(rf_ctx *c, int n, u64 modulus, u64 a, u64 b, int k, const u64 *v, int w,
             const u64 *Qimg, const u64 *sqrtimg, const u64 *chkimg, const u64 *piimg, u64 trmask)
{
    memset(c, 0, sizeof *c);
    c->n = n; c->k = k; c->w = w;
    c->mask = (n == 64) ? ~0ULL : ((1ULL << n) - 1);
    c->red = modulus & c->mask;
    c->a = a; c->b = b; c->trmask = trmask;
    lm_build(&c->Q, Qimg, n);
    lm_build(&c->sqrt, sqrtimg, n);
    lm_build(&c->chk, chkimg, n);
    lm_build(&c->pi, piimg, n);
    c->sqrtb = lm(&c->sqrt, b);
    for (int i = 0; i < k; i++) {
        c->v[i] = v[i];
        c->wv[i] = lm(&c->Q, v[i]) ^ v[i];
        c->trv |= (u64)ftr(c, v[i]) << i;
    }
}

size_t rf_ctx_size(void) { return sizeof(rf_ctx); }

/* arithmetic exports, for the cross-check against the Python field */
u64 rf_mul(const rf_ctx *c, u64 x, u64 y) { return fmul(c, x, y); }
u64 rf_inv(const rf_ctx *c, u64 x) { return finv(c, x); }
void rf_add(const rf_ctx *c, const u64 *P, const u64 *Q, u64 *R)
{
    pt p = {P[0], P[1], (int)P[2]}, q = {Q[0], Q[1], (int)Q[2]}, r = padd(c, p, q);
    R[0] = r.x; R[1] = r.y; R[2] = (u64)r.inf;
}

/* Factor base: every point with abscissa in V (both signs; the 2-torsion point once). Returns the count. */
long rf_fb(const rf_ctx *c, u64 *xs, u64 *ys, long cap)
{
    long cnt = 0;
    u64 x = 0;
    for (u64 g = 0; g < (1ULL << c->k); g++) {
        if (g) x ^= c->v[__builtin_ctzll(g)];           /* Gray code walk over V */
        if (x == 0) {
            if (cnt < cap) { xs[cnt] = 0; ys[cnt] = c->sqrtb; }
            cnt++;
            continue;
        }
        u64 ix = finv(c, x);
        u64 cc = x ^ c->a ^ fmul(c, c->b, fsq(c, ix));
        if (ftr(c, cc)) continue;
        u64 y = fmul(c, x, lm(&c->Q, cc));
        if (cnt + 1 < cap) { xs[cnt] = x; ys[cnt] = y; xs[cnt + 1] = x; ys[cnt + 1] = y ^ x; }
        cnt += 2;
    }
    return cnt;
}

#define QCAP 256
#define QPTS 6
/* a solution of the linear system, waiting for the quadratic x^2 + e1 x + e2 (inversions of e1 are shared) */
typedef struct { u64 Tx, ym, e1, e2, lam; int nchosen; pt chosen[QPTS]; } cand;

typedef struct {
    const rf_ctx *c;
    int m;
    const u64 *fx, *fy;
    long nfb;
    u64 *out;
    long cap, cnt;
    rf_stats *st;
    pt chosen[16];
    cand q[QCAP];
    int nq;
} run;

static void emit(run *r, int nchosen, pt P1, pt P2)
{
    if (r->cnt < r->cap) {
        u64 *o = r->out + (size_t)r->cnt * 2 * r->m;
        for (int i = 0; i < nchosen; i++) { *o++ = r->chosen[i].x; *o++ = r->chosen[i].y; }
        *o++ = P1.x; *o++ = P1.y;
        if (nchosen + 2 == r->m) { *o++ = P2.x; *o++ = P2.y; }
    }
    r->cnt++;
}

/* Two-summand step by enumeration: every factor-base P1 with x(T - P1) in V (index order P1 <= T - P1 not
 * imposed; the caller deduplicates).  Used as the baseline leaf and for the target with x(T) = 0. */
static void d2_enum(run *r, int nchosen, pt T)
{
    const rf_ctx *c = r->c;
    for (long i = 0; i < r->nfb; i++) {
        pt P = {r->fx[i], r->fy[i], 0};
        pt Q = padd(c, T, pneg(P));
        if (Q.inf) continue;
        if (lm(&c->chk, Q.x)) continue;
        emit(r, nchosen, P, Q);
    }
}

static void batch_inv(const rf_ctx *c, const u64 *in, u64 *out, int cnt);

static void flush(run *r)
{
    const rf_ctx *c = r->c;
    rf_stats *st = r->st;
    u64 e1[QCAP], ie[QCAP];
    pt save[QPTS];
    memcpy(save, r->chosen, sizeof save);
    for (int i = 0; i < r->nq; i++) e1[i] = r->q[i].e1;
    batch_inv(c, e1, ie, r->nq);
    for (int i = 0; i < r->nq; i++) {
        const cand *q = &r->q[i];
        u64 x1, x2;
        if (q->e1 == 0) {
            x1 = x2 = lm(&c->sqrt, q->e2);
        } else {
            u64 d = fmul(c, q->e2, fsq(c, ie[i]));
            if (ftr(c, d)) continue;
            x1 = fmul(c, q->e1, lm(&c->Q, d));
            x2 = x1 ^ q->e1;
        }
        st->quad_split++;
        if (lm(&c->chk, x1)) continue;                       /* x2 = x1 + e1 is then in V as well */
        u64 mu = fmul(c, q->lam, q->Tx) ^ q->ym;
        /* P1, P2 and -T are the three points of the line, so P1 + P2 = T */
        pt P1 = {x1, fmul(c, q->lam, x1) ^ mu, 0}, P2 = {x2, fmul(c, q->lam, x2) ^ mu, 0};
        st->found++;
        for (int j = 0; j < q->nchosen; j++) r->chosen[j] = q->chosen[j];
        emit(r, q->nchosen, P1, P2);
    }
    memcpy(r->chosen, save, sizeof save);
    r->nq = 0;
}

/* Two-summand step by linear algebra.  P1 + P2 = T  <=>  P1, P2, -T on a line y = lam (x + xT) + yT + xT.
 * Then  e1 = x1 + x2 = lam^2 + lam + a + xT   and   e2 = x1 x2 = xT lam^2 + ((yT + xT)^2 + b) / xT,
 * both F_2-affine in lam.  With e1 = sum t_i v_i in V:  lam = Q(e1 + c1) + eps  (needs Tr(e1 + c1) = 0), and
 * e2 = xT W(e1) + eps xT + xT W(c1) + c2  must lie in V2 = span(V V): a linear system in (t, eps), k + 1 unknowns,
 * n - dim V2 + 1 equations.  Each solution gives x^2 + e1 x + e2, whose roots must lie in V. */
static void d2_core(run *r, int nchosen, pt T, u64 ixT)
{
    const rf_ctx *c = r->c;
    const int k = c->k;
    rf_stats *st = r->st;
    st->targets++;
    u64 ym = T.y ^ T.x;                         /* y(-T) */
    u64 c1 = c->a ^ T.x;
    u64 c2 = fmul(c, fsq(c, ym) ^ c->b, ixT);
    const u64 trbit = 1ULL << (c->w - 1);
    u64 rhs = lm(&c->pi, fmul(c, T.x, lm(&c->Q, c1) ^ c1) ^ c2) | (ftr(c, c1) ? trbit : 0);

    /* xor basis of the columns, with the combination of unknowns that produced each element */
    u64 bvec[64], bcomb[64], kern[65];
    int piv[64], nb = 0, nk = 0;
    for (int i = 0; i <= k; i++) {
        u64 col = (i < k) ? (lm(&c->pi, fmul(c, T.x, c->wv[i])) | (((c->trv >> i) & 1) ? trbit : 0))
                          : lm(&c->pi, T.x);
        u64 comb = 1ULL << i;
        for (int j = 0; j < nb; j++) {                      /* branch-free: the bits are unpredictable */
            u64 mk = -((col >> piv[j]) & 1);
            col ^= bvec[j] & mk; comb ^= bcomb[j] & mk;
        }
        if (col) { piv[nb] = 63 - __builtin_clzll(col); bvec[nb] = col; bcomb[nb] = comb; nb++; }
        else kern[nk++] = comb;
    }
    u64 part = 0;
    for (int j = 0; j < nb; j++) {
        u64 mk = -((rhs >> piv[j]) & 1);
        rhs ^= bvec[j] & mk; part ^= bcomb[j] & mk;
    }
    if (rhs) return;
    st->lin_solvable++;
    if ((u64)nk > st->max_kernel_dim) st->max_kernel_dim = (u64)nk;

    u64 sol = part;
    for (u64 g = 0; g < (1ULL << nk); g++) {
        if (g) sol ^= kern[__builtin_ctzll(g)];
        st->lin_candidates++;
        u64 e1 = 0, tm = sol & ((1ULL << k) - 1), eps = (sol >> k) & 1;
        while (tm) { e1 ^= c->v[__builtin_ctzll(tm)]; tm &= tm - 1; }
        cand *q = &r->q[r->nq++];
        q->Tx = T.x; q->ym = ym; q->e1 = e1;
        q->lam = lm(&c->Q, e1 ^ c1) ^ eps;
        q->e2 = fmul(c, T.x, q->lam ^ e1 ^ c1) ^ c2;         /* lam^2 = lam + e1 + c1 */
        q->nchosen = nchosen;
        for (int i = 0; i < nchosen; i++) q->chosen[i] = r->chosen[i];
        if (r->nq == QCAP) flush(r);
    }
}

static void d2_linear(run *r, int nchosen, pt T)
{
    if (T.x == 0) { r->st->targets++; r->st->degenerate++; d2_enum(r, nchosen, T); return; }
    d2_core(r, nchosen, T, finv(r->c, T.x));
}

/* out[i] = 1 / in[i] for in[i] != 0 (zeros are skipped and left 0): one inversion, 3 multiplications each */
static void batch_inv(const rf_ctx *c, const u64 *in, u64 *out, int cnt)
{
    u64 acc = 1;
    for (int i = 0; i < cnt; i++) { out[i] = acc; if (in[i]) acc = fmul(c, acc, in[i]); }
    acc = finv(c, acc);
    for (int i = cnt - 1; i >= 0; i--) {
        if (!in[i]) { out[i] = 0; continue; }
        out[i] = fmul(c, out[i], acc);
        acc = fmul(c, acc, in[i]);
    }
}

#define BLK 256
/* The innermost enumerated point, in blocks with shared inversions: T_i = S - P_i, then the leaf on T_i. */
static void last_level(run *r, int depth, long start, long stop, pt S, int linear)
{
    const rf_ctx *c = r->c;
    u64 d[BLK], id[BLK], tx[BLK], itx[BLK];
    pt T[BLK];
    for (long base = start; base < stop; base += BLK) {
        int cnt = (int)(stop - base < BLK ? stop - base : BLK);
        for (int j = 0; j < cnt; j++) d[j] = S.inf ? 0 : (S.x ^ r->fx[base + j]);
        batch_inv(c, d, id, cnt);
        for (int j = 0; j < cnt; j++) {
            u64 px = r->fx[base + j], py = r->fy[base + j];
            if (d[j] == 0) {                                  /* S = O, or S = +-P: generic addition */
                pt P = {px, py, 0};
                T[j] = padd(c, S, pneg(P));
            } else {
                u64 lam = fmul(c, S.y ^ py ^ px, id[j]);
                T[j].x = fsq(c, lam) ^ lam ^ d[j] ^ c->a;
                T[j].y = fmul(c, lam, S.x ^ T[j].x) ^ T[j].x ^ S.y;
                T[j].inf = 0;
            }
            tx[j] = T[j].inf ? 0 : T[j].x;
        }
        if (linear) batch_inv(c, tx, itx, cnt);
        for (int j = 0; j < cnt; j++) {
            if (T[j].inf) continue;
            pt P = {r->fx[base + j], r->fy[base + j], 0};
            r->chosen[depth] = P;
            if (!linear) {
                if (lm(&c->chk, T[j].x) == 0) { pt none = {0, 0, 0}; emit(r, depth + 1, T[j], none); }
            } else if (tx[j] == 0) d2_linear(r, depth + 1, T[j]);
            else d2_core(r, depth + 1, T[j], itx[j]);
        }
    }
}

static void rec(run *r, int depth, int nguess, long start, pt S, int linear)
{
    const rf_ctx *c = r->c;
    if (depth == nguess - 1) { last_level(r, depth, start, r->nfb, S, linear); return; }
    if (depth == nguess) {
        if (S.inf) return;
        if (linear) d2_linear(r, depth, S);
        else if (lm(&c->chk, S.x) == 0) { pt none = {0, 0, 0}; emit(r, depth, S, none); }
        return;
    }
    for (long i = start; i < r->nfb; i++) {
        pt P = {r->fx[i], r->fy[i], 0};
        r->chosen[depth] = P;
        rec(r, depth + 1, nguess, i, padd(c, S, pneg(P)), linear);
    }
}

/* Decompositions R = P_1 + ... + P_m with every x(P_i) in V.  Rows of 2m words (x, y per point); the same
 * multiset can appear more than once.  Outer index range [lo, hi) of the first enumerated point allows the
 * caller to split the work.  Returns the number of rows found (rows beyond cap are counted, not stored). */
long rf_lrf(const rf_ctx *c, int m, const u64 *R, const u64 *fx, const u64 *fy, long nfb,
            long lo, long hi, u64 *out, long cap, rf_stats *st)
{
    static __thread run r;
    r.c = c; r.m = m; r.fx = fx; r.fy = fy; r.nfb = nfb; r.out = out; r.cap = cap; r.cnt = 0; r.st = st; r.nq = 0;
    if (m - 2 > QPTS) return -2;
    pt S = {R[0], R[1], 0};
    memset(st, 0, sizeof *st);
    if (hi > nfb) hi = nfb;
    if (c->w - 1 <= 1) {                  /* dim V2 >= n - 1: the linear system says nothing, enumerate */
        st->degenerate = 1;
        if (m == 2) { last_level(&r, 0, lo, hi, S, 0); return r.cnt; }
        for (long i = lo; i < hi; i++) {
            pt P = {fx[i], fy[i], 0};
            r.chosen[0] = P;
            rec(&r, 1, m - 1, i, padd(c, S, pneg(P)), 0);
        }
        return r.cnt;
    }
    if (m == 2) { d2_linear(&r, 0, S); flush(&r); return r.cnt; }
    if (m == 3) { last_level(&r, 0, lo, hi, S, 1); flush(&r); return r.cnt; }
    for (long i = lo; i < hi; i++) {
        pt P = {fx[i], fy[i], 0};
        r.chosen[0] = P;
        rec(&r, 1, m - 2, i, padd(c, S, pneg(P)), 1);
    }
    flush(&r);
    return r.cnt;
}

long rf_brute(const rf_ctx *c, int m, const u64 *R, const u64 *fx, const u64 *fy, long nfb,
              long lo, long hi, u64 *out, long cap)
{
    rf_stats st;
    static __thread run r;
    r.c = c; r.m = m; r.fx = fx; r.fy = fy; r.nfb = nfb; r.out = out; r.cap = cap; r.cnt = 0; r.st = &st; r.nq = 0;
    pt S = {R[0], R[1], 0};
    memset(&st, 0, sizeof st);
    if (hi > nfb) hi = nfb;
    if (m == 2) { last_level(&r, 0, lo, hi, S, 0); return r.cnt; }
    for (long i = lo; i < hi; i++) {
        pt P = {fx[i], fy[i], 0};
        r.chosen[0] = P;
        rec(&r, 1, m - 1, i, padd(c, S, pneg(P)), 0);
    }
    return r.cnt;
}

/* m = 4 meet in the middle: sorted table of (x, y, i, j) for P_i + P_j, i <= j; lookups of R - P_k - P_l. */
typedef struct { u64 x, y; uint32_t i, j; } pair;
static int cmp_pair(const void *A, const void *B)
{
    const pair *p = A, *q = B;
    if (p->x != q->x) return p->x < q->x ? -1 : 1;
    if (p->y != q->y) return p->y < q->y ? -1 : 1;
    return 0;
}

long rf_mitm4(const rf_ctx *c, const u64 *R, const u64 *fx, const u64 *fy, long nfb, u64 *out, long cap,
              u64 *table_bytes)
{
    size_t np = (size_t)nfb * (nfb + 1) / 2, n = 0;
    pair *tab = malloc(np * sizeof *tab);
    if (!tab) return -1;
    *table_bytes = np * sizeof *tab;
    for (long i = 0; i < nfb; i++)
        for (long j = i; j < nfb; j++) {
            pt P = {fx[i], fy[i], 0}, Q = {fx[j], fy[j], 0}, s = padd(c, P, Q);
            if (s.inf) continue;
            tab[n].x = s.x; tab[n].y = s.y; tab[n].i = (uint32_t)i; tab[n].j = (uint32_t)j; n++;
        }
    qsort(tab, n, sizeof *tab, cmp_pair);
    long cnt = 0;
    pt S = {R[0], R[1], 0};
    for (long i = 0; i < nfb; i++)
        for (long j = i; j < nfb; j++) {
            pt P = {fx[i], fy[i], 0}, Q = {fx[j], fy[j], 0};
            pt t = padd(c, padd(c, S, pneg(P)), pneg(Q));
            if (t.inf) continue;
            size_t a = 0, b = n;
            while (a < b) {
                size_t mid = (a + b) / 2;
                if (tab[mid].x < t.x || (tab[mid].x == t.x && tab[mid].y < t.y)) a = mid + 1; else b = mid;
            }
            for (; a < n && tab[a].x == t.x && tab[a].y == t.y; a++) {
                if (cnt < cap) {
                    u64 *o = out + (size_t)cnt * 8;
                    o[0] = P.x; o[1] = P.y; o[2] = Q.x; o[3] = Q.y;
                    o[4] = fx[tab[a].i]; o[5] = fy[tab[a].i]; o[6] = fx[tab[a].j]; o[7] = fy[tab[a].j];
                }
                cnt++;
            }
        }
    free(tab);
    return cnt;
}
