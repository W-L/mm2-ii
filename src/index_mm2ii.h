#ifndef MM2_II_INDEX_MM2II_H
#define MM2_II_INDEX_MM2II_H

#include "minimap.h"

// copy from static def in mm2 index.c
typedef struct mm_idx_bucket_s {
    mm128_v a;   // (minimizer, position) array
    int32_t n;   // size of the _p_ array
    uint64_t *p; // position array for minimizers appearing >1 times
    void *h;     // hash table indexing _p_ and minimizers appearing once
} mm_idx_bucket_t;


mm_idx_t *mm_idx_init(int w, int k, int b, int flag);   // pulling from mm2 index.c

mm_idx_t *init_index(int w, int k);                     // wrapper for mm_idx_init in mm2 index.c

void mm_idx_add(mm_idx_t *mi, int n, const mm128_t *a); // copy defd in index_mm2ii.c

void mm_idx_dump(FILE *fp, const mm_idx_t *mi);         // pulling from mm2 index.c

void mm_idx_post(mm_idx_t *mi, int n_threads);

void mm_idx_destroy(mm_idx_t *mi);


#endif //MM2_II_INDEX_MM2II_H
