/* In-place GF(2) echelon (rank + pivot columns) of a bit-packed matrix.
   Blocked elimination, 8 pivots at a time (method of the four Russians on the pivot block), OpenMP over rows.
   Convention as csd/_elim.c: columns scanned left to right, bit c of a row = word c/64, bit c%64.
   The pivot COLUMNS returned equal those of the plain row echelon (first column order), tested against it.
   The matrix is destroyed (memory = rows*words*8 plus a 256-row table). */
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>

#define K 8
static inline int getbit(const uint64_t *row, long c) { return (row[c >> 6] >> (c & 63)) & 1; }

long gf2_rank_inplace(uint64_t *M, long rows, long words, long cols, long *piv, int nthreads)
{
    long rank = 0;
    omp_set_num_threads(nthreads);
    uint64_t *T = (uint64_t *)malloc((size_t)(1 << K) * words * 8);
    if (!T) return -1;
    long c = 0;
    while (c < cols && rank < rows) {
        long pc[K]; int cnt = 0;
        long w0 = c >> 6;
        for (; c < cols && cnt < K && rank + cnt < rows; c++) {
            long found = -1;
            for (long r = rank + cnt; r < rows; r++) {
                uint64_t *row = M + r * words;
                int b = getbit(row, c);
                int eff[K];
                for (int j = 0; j < cnt; j++) {
                    int e = getbit(row, pc[j]);
                    for (int i = 0; i < j; i++) if (eff[i]) e ^= getbit(M + (rank + i) * words, pc[j]);
                    eff[j] = e;
                    if (e) b ^= getbit(M + (rank + j) * words, c);
                }
                if (b) { found = r; break; }
            }
            if (found < 0) continue;
            long dst = rank + cnt;
            uint64_t *pr = M + found * words;
            /* reduce the chosen row by earlier block pivots (in order) */
            for (int j = 0; j < cnt; j++)
                if (getbit(pr, pc[j])) { uint64_t *q = M + (rank + j) * words; for (long k = w0; k < words; k++) pr[k] ^= q[k]; }
            if (found != dst) {
                uint64_t *b2 = M + dst * words;
                for (long k = w0; k < words; k++) { uint64_t t = pr[k]; pr[k] = b2[k]; b2[k] = t; }
            }
            pc[cnt++] = c;
            piv[rank + cnt - 1] = c;
        }
        if (cnt == 0) break;
        /* make the block pivot rows reduced w.r.t. each other */
        for (int j = cnt - 1; j >= 0; j--) {
            uint64_t *q = M + (rank + j) * words;
            for (int i = 0; i < j; i++) {
                uint64_t *p = M + (rank + i) * words;
                if (getbit(p, pc[j])) for (long k = w0; k < words; k++) p[k] ^= q[k];
            }
        }
        /* table of all combinations */
        memset(T, 0, (size_t)words * 8);
        for (int idx = 1; idx < (1 << cnt); idx++) {
            int low = __builtin_ctz(idx);
            uint64_t *t = T + (size_t)idx * words;
            const uint64_t *prev = T + (size_t)(idx & (idx - 1)) * words;
            const uint64_t *q = M + (rank + low) * words;
            for (long k = w0; k < words; k++) t[k] = prev[k] ^ q[k];
        }
        long r0 = rank + cnt;
        #pragma omp parallel for schedule(static, 64)
        for (long i = r0; i < rows; i++) {
            uint64_t *ri = M + i * words;
            int idx = 0;
            for (int j = 0; j < cnt; j++) idx |= getbit(ri, pc[j]) << j;
            if (idx) {
                const uint64_t *t = T + (size_t)idx * words;
                for (long k = w0; k < words; k++) ri[k] ^= t[k];
            }
        }
        rank += cnt;
    }
    free(T);
    return rank;
}
