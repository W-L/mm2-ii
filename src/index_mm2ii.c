//#include <stdio.h>

#include "index_mm2ii.h"

#include "minimap.h"
#include "kvec.h"
#include "khash.h"
#include "kthread.h"
#include "mmpriv.h"

#define idx_hash(a) ((a)>>1)
#define idx_eq(a, b) ((a)>>1 == (b)>>1)
KHASH_INIT(idx, uint64_t, uint64_t, 1, idx_hash, idx_eq)
typedef khash_t(idx) idxhash_t;
KHASH_MAP_INIT_STR(str, uint32_t)


// wrapper for mm_idx_init
mm_idx_t *init_index(int w, int k)
{
    mm_idx_t *mi;
    int bucket_bits = 14;  // these are default from index.c mm_idx_str TODO check if this makes diff
    int flag = 0;
    // definition in mm2's index.c
    mi = mm_idx_init(w, k, bucket_bits, flag);
    return mi;
}



// Function definitions copied from mm2's index.c
// These are declared static in mm2's source code, so can not be included directly

void mm_idx_add(mm_idx_t *mi, int n, const mm128_t *a)
{
    int i, mask = (1<<mi->b) - 1;
    for (i = 0; i < n; ++i) {
        mm128_v *p = &mi->B[a[i].x>>8&mask].a;
        kv_push(mm128_t, 0, *p, a[i]);
    }
}


void worker_post(void *g, long i, int tid)  // unused parameter for kt_for
{
    int n, n_keys;
    size_t j, start_a, start_p;
    idxhash_t *h;
    mm_idx_t *mi = (mm_idx_t*)g;
    mm_idx_bucket_t *b = &mi->B[i];
    if (b->a.n == 0) return;

    // sort by minimizer
    radix_sort_128x(b->a.a, b->a.a + b->a.n);

    // count and preallocate
    for (j = 1, n = 1, n_keys = 0, b->n = 0; j <= b->a.n; ++j) {
        if (j == b->a.n || b->a.a[j].x>>8 != b->a.a[j-1].x>>8) {
            ++n_keys;
            if (n > 1) b->n += n;
            n = 1;
        } else ++n;
    }
    h = kh_init(idx);
    kh_resize(idx, h, n_keys);
    b->p = (uint64_t*)calloc(b->n, 8);

    // create the hash table
    for (j = 1, n = 1, start_a = start_p = 0; j <= b->a.n; ++j) {
        if (j == b->a.n || b->a.a[j].x>>8 != b->a.a[j-1].x>>8) {
            khint_t itr;
            int absent;
            mm128_t *p = &b->a.a[j-1];
            itr = kh_put(idx, h, p->x>>8>>mi->b<<1, &absent);
            assert(absent && j == start_a + n);
            if (n == 1) {
                kh_key(h, itr) |= 1;
                kh_val(h, itr) = p->y;
            } else {
                int k;
                for (k = 0; k < n; ++k)
                    b->p[start_p + k] = b->a.a[start_a + k].y;
                radix_sort_64(&b->p[start_p], &b->p[start_p + n]); // sort by position; needed as in-place radix_sort_128x() is not stable
                kh_val(h, itr) = (uint64_t)start_p<<32 | n;
                start_p += n;
            }
            start_a = j, n = 1;
        } else ++n;
    }
    b->h = h;
    assert(b->n == (int32_t)start_p);

    // deallocate and clear b->a
    kfree(0, b->a.a);
    b->a.n = b->a.m = 0, b->a.a = 0;
}


void mm_idx_post(mm_idx_t *mi, int n_threads)
{
    // worker_post(mi, 0, 0);  // debugging: call on single bucket
    kt_for(n_threads, worker_post, mi, 1<<mi->b);
}
