/* Forward Gaussian elimination kernels for csd/elim.py.
 *
 * Both kernels reduce a dense matrix to row echelon form in place, scanning the
 * columns left to right, and write the pivot column of each pivot row. The
 * caller orders the columns so that the pivots it cares about (low-degree
 * monomials, the constant) come last: the number of pivots in a final block of
 * columns is then the dimension of the row-space vectors supported on it.
 *
 * gf2_echelon: rows are bit-packed, `words` uint64 per row, bit c of the row is
 *              word c / 64, bit c % 64.
 * fp_echelon:  entries are uint32 residues mod p, p < 2^31.
 * Each returns the rank and counts the row operations it performed in *ops
 * (word XORs for GF(2), multiply-adds for F_p): counted, not timed.
 */
#include <stdint.h>
#include <string.h>
#include <stdlib.h>

long gf2_echelon(uint64_t *M, long rows, long words, long cols, long *piv, double *ops)
{
    long rank = 0;
    double cnt = 0;
    for (long c = 0; c < cols && rank < rows; c++) {
        long w = c >> 6;
        uint64_t bit = (uint64_t)1 << (c & 63);
        long r = rank;
        while (r < rows && !(M[r * words + w] & bit)) r++;
        if (r == rows) continue;
        if (r != rank) {
            uint64_t *a = M + r * words, *b = M + rank * words;
            for (long k = w; k < words; k++) { uint64_t t = a[k]; a[k] = b[k]; b[k] = t; }
        }
        uint64_t *pr = M + rank * words;
        for (long i = rank + 1; i < rows; i++) {
            uint64_t *ri = M + i * words;
            if (ri[w] & bit) {
                for (long k = w; k < words; k++) ri[k] ^= pr[k];
                cnt += (double)(words - w);
            }
        }
        piv[rank++] = c;
    }
    *ops = cnt;
    return rank;
}

static uint32_t inv_mod(uint32_t a, uint32_t p)
{
    int64_t t = 0, nt = 1, r = p, nr = a;
    while (nr) {
        int64_t q = r / nr, x;
        x = t - q * nt; t = nt; nt = x;
        x = r - q * nr; r = nr; nr = x;
    }
    if (t < 0) t += p;
    return (uint32_t)t;
}

static inline uint32_t red(uint64_t x, uint64_t p, double invp)
{
    /* x mod p for x < 2^62 by a float reciprocal, corrected exactly */
    int64_t q = (int64_t)((double)x * invp);
    int64_t r = (int64_t)x - q * (int64_t)p;
    while (r < 0) r += p;
    while (r >= (int64_t)p) r -= p;
    return (uint32_t)r;
}

long fp_echelon(uint32_t *M, long rows, long cols, uint32_t p, long *piv, double *ops)
{
    long rank = 0;
    double cnt = 0;
    double invp = 1.0 / (double)p;
    long *nzk = (long *)malloc(sizeof(long) * (cols > 0 ? cols : 1));
    for (long c = 0; c < cols && rank < rows; c++) {
        long r = rank;
        while (r < rows && M[r * cols + c] == 0) r++;
        if (r == rows) continue;
        if (r != rank) {
            uint32_t *a = M + r * cols, *b = M + rank * cols;
            for (long k = c; k < cols; k++) { uint32_t t = a[k]; a[k] = b[k]; b[k] = t; }
        }
        uint32_t *pr = M + rank * cols;
        uint64_t iv = inv_mod(pr[c], p);
        long nnz = 0;
        for (long k = c; k < cols; k++) {
            if (pr[k]) {
                pr[k] = red((uint64_t)pr[k] * iv, p, invp);
                nzk[nnz++] = k;
            }
        }
        for (long i = rank + 1; i < rows; i++) {
            uint32_t *ri = M + i * cols;
            uint64_t f = ri[c];
            if (f) {
                uint64_t g = p - f;
                for (long t = 0; t < nnz; t++) {
                    long k = nzk[t];
                    ri[k] = red((uint64_t)ri[k] + g * pr[k], p, invp);
                }
                cnt += (double)(cols - c);  /* counted as a dense row operation */
            }
        }
        piv[rank++] = c;
    }
    free(nzk);
    *ops = cnt;
    return rank;
}
