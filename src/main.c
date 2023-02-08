#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>
#include <time.h>

// external dependencies
#include <zlib.h>

// minimap header dependencies
#include "minimap.h"
#include "khash.h"
#include "kseq.h"
#include "mmpriv.h"     // radix_sort, mm_seq4..

//#include "index_mm2ii.h"

// inits for kseq & khash
KSEQ_INIT(gzFile, gzread)

KHASH_SET_INIT_INT(int_set)
typedef khash_t(int_set) intset_t;

KHASH_MAP_INIT_STR(str, uint32_t)
typedef khash_t(str) strdict_t;


// magic identifier for mmis
#define MM2_VERSION     "2.24-r1150"
#define MM2_II_VERSION  "0.0.3"
#define WINDOW_SIZE     5
#define KMER_SIZE       15
//#define NTHREADS        6



// main struct to hold minimizers & sequence name
typedef struct{
    char n[64];
    mm128_v v;
} Sketch;


typedef struct{
    uint32_t x[5];
} idx_info;


typedef struct{
    int *a;
    int l;
} int_arr;


typedef struct{
    int n_seq, sum_len, header_len;
} seq_counter;


// REMINDER
//typedef struct mm_idx_bucket_s {
//    mm128_v a;   // (minimizer, position) array
//    int32_t n;   // size of the _p_ array
//    uint64_t *p; // position array for minimizers appearing >1 times
//    void *h;     // hash table indexing _p_ and minimizers appearing once
//} mm_idx_bucket_t;


typedef struct {
    uint32_t s;
    int32_t n;   // size of the _p_ array
    uint64_t *p; // position array for minimizers appearing >1 times
    uint64_t *keys;
    uint64_t *vals;
} idx_bucket;



int in_int_set(intset_t *ih, const uint32_t query)
{
    khint_t q;
    int in_set;
    q = kh_get(int_set, ih, query);
    in_set = (q != kh_end(ih));
    if (in_set) return 1;
    return 0;
}



void dump_sequence_info(FILE *fp_in, FILE *fp_out, intset_t *idx_del, uint32_t n_seq) {
    for (int i = 0; i < n_seq; ++i) {
        uint8_t l;
        fread(&l, 1, 1, fp_in);
        fwrite(&l, 1, 1, fp_out);
        if (l) {
            char name[l];
            fread(&name, 1, l, fp_in);
            fwrite(&name, 1, l, fp_out);
        }
        uint32_t len;
        fread(&len, 4, 1, fp_in);
        // length of deleted seqs set to 0
        if (in_int_set(idx_del, i)){
            len = 0;
        }
        fwrite(&len, 4, 1, fp_out);
    }
}



void dump_basics(FILE *fp_in, FILE *fp_out, idx_info *mm_info)
{
    char magic[4];
    fread(magic, 1, 4, fp_in);
    fwrite(magic, 1, 4, fp_out);
    // w, k, b, n_seq, flag
    fread(mm_info->x, 4, 5, fp_in);
    fwrite(mm_info->x, 4, 5, fp_out);
}



uint64_t *reduce_to_singleton(idx_bucket *b, int ppos, int mm_idx)
{
    uint64_t mm[2];
    mm[0] = b->keys[mm_idx], mm[1] = b->vals[mm_idx];
    // new minimizer
    uint64_t *nn = malloc(sizeof(uint64_t) * 2);
    nn[0] = mm[0]|1, nn[1] = b->p[ppos];
    return nn;
}



uint64_t *change_modified_multiton(idx_bucket *b, int pshift, int newsize, int mm_idx)
{
    uint64_t mm[2];
    mm[0] = b->keys[mm_idx], mm[1] = b->vals[mm_idx];
    // new minimizer value with mod to pshift and copynumber
    uint32_t ppos = mm[1] >> 32;
    ppos = ppos - pshift;
    uint64_t *nn = malloc(sizeof(uint64_t) * 2);
    nn[0] = mm[0], nn[1] = (uint64_t) ppos<<32 | newsize;
    return nn;
}



int array_sum(const int *arr, int size){
    int sum = 0;
    for (int i = 0; i < size; i++){
        sum += arr[i];
    }
    return sum;
}



int *array_cumsum(const int *arr, int size){
    int *cumsum = malloc(sizeof(int) * size);
    int sum = 0;
    for (int i = 0; i < size; i++){
        sum += arr[i];
        cumsum[i] = sum;
    }
    return cumsum;
}



void modify_parray(FILE *fp, idx_bucket *b, int *pdel)
{
    int n_del_p = array_sum(pdel, b->n);
    int32_t new_p_size = 0;
    new_p_size = b->n - n_del_p;
    assert(new_p_size != 1); // reduction should delete entire p array
    assert(new_p_size >= 0); // can't be negative

    if (n_del_p == 0) {   // no p to delete
        fwrite(&b->n, 4, 1, fp);
        fwrite(b->p, 8, b->n, fp);
        return;
    }
    // elements in p to delete
    // no more multitons left, p array is now empty
    if (new_p_size == 0) {
        fwrite(&new_p_size, 4, 1, fp);
        fwrite(b->p, 8, new_p_size, fp);
        return;
    }
    // some elements of p array remain
    fwrite(&new_p_size, 4, 1, fp);
    for (int m = 0; m < b->n; m++) {
        if (pdel[m]) {continue;}
        fwrite(&b->p[m], 8, 1, fp);
    }
}



int is_multiton(uint64_t mm)
{
    if ((mm & 1) != 1){
        return 1;
    } else {
        return 0;
    }
}



int *process_bucket(idx_bucket *b, intset_t *idx_del, uint32_t *n_del_s, uint32_t *n_del_m, uint32_t *n_red, uint32_t *n_mod){
    int m = b->n;

    // helper arrays for multitons mods
    int cn_arr[m];
    int idx_arr[m];
    int newsize_arr[m];
    int *pdel = malloc(sizeof(int) * m);

    for (int i = 0; i < b->s; i++) {
        uint64_t mm[2];
        mm[0] = b->keys[i], mm[1] = b->vals[i];

        if (is_multiton(mm[0])) {               // multiton
            // get copy number
            uint32_t copyn;
            copyn = (uint32_t) mm[1];
            uint32_t ppos = mm[1] >> 32;
            int n_del_units = 0;

            // fill the respective positions in karr and carr
            for (int l = 0; l < copyn; l++) {
                cn_arr[ppos + l] = copyn;
                idx_arr[ppos + l] = i;

                uint64_t seq_id = b->p[ppos + l] >> 32;
                // check if source sequence is deleted
                if (in_int_set(idx_del, seq_id)) {
                    pdel[ppos + l] = 1;
                    n_del_units++;
                } else {
                    pdel[ppos + l] = 0;
                }
            }

            // check if still multiton here, otherwise mark all p units for deletion
            if (copyn - n_del_units == 0){
                for (int r = 0; r < copyn; r++) {
                    pdel[ppos + r] = 1;
                }
            }

            for (int r = 0; r < copyn; r++) {
                newsize_arr[ppos + r] = copyn - n_del_units;
            }

        } else {                                    // singleton
            // check if source sequence is deleted
            if (in_int_set(idx_del, mm[1] >> 32)) {
                // delete mm
                (*n_del_s)++;
                b->keys[i] = 0;
                b->vals[i] = 0;
            }

        }
    }


    // effect multiton info arrays
    for (int l = 0; l < m;) {
        int cn = cn_arr[l];
        int idx = idx_arr[l];
        int newsize = newsize_arr[l];
        assert(newsize >= 0);

        // no change
        if (cn == newsize){
            l += cn;
            continue;}

        // wipe-out
        if (newsize == 0){
            (*n_del_m)++;
            b->keys[idx] = 0;
            b->vals[idx] = 0;
            l += cn;
            continue;
        }

        // reduce to singleton
        if (newsize == 1){
            // search for the copy in p that comes from a non-deleted sequence
            int done = 0;
            uint64_t *nn;
            for (int i = 0; i < cn; i++) {
                if (pdel[l + i] == 0) {
                    assert(!done);
                    nn = reduce_to_singleton(b, l + i, idx);
                    b->keys[idx] = nn[0];
                    b->vals[idx] = nn[1];
                    pdel[l + i] = 1;
                    done = 1;
                }
            }
            assert(done);
            (*n_red)++;
            l += cn;
            continue;
        }

        // remain multiton
        if (newsize > 1){
            // get pshift
            int *pdel_cumsum;
            pdel_cumsum = array_cumsum(pdel, m);
            int pshift = pdel_cumsum[l];
            uint64_t *nn = change_modified_multiton(b, pshift, newsize, idx);
            (*n_mod)++;
            b->keys[idx] = nn[0];
            b->vals[idx] = nn[1];
            l += cn;
            continue;
        }
    }
    return pdel;
}



void dump_minimizers(FILE *fp_out, idx_bucket *b, const uint32_t *n_del_s, const uint32_t *n_del_m)
{
    // empty bucket
    if (b->s == 0){
        fwrite(&b->s, 4, 1, fp_out);
        return;
    }

    // non-empty bucket
    int del = 0;
    uint32_t size = b->s - *n_del_s - *n_del_m;
    fwrite(&size, 4, 1, fp_out);
    for (int i = 0; i < b->s; i++) {
        if ((b->keys[i] == 0) && (b->vals[i] == 0)){
            del++;
            continue;
        }
        fwrite(&b->keys[i], 8, 1, fp_out);
        fwrite(&b->vals[i], 8, 1, fp_out);
    }
    assert(del == (*n_del_s + *n_del_m));
}



seq_counter count_sequences(const char * fname){
    gzFile file = gzopen(fname, "r");
    kseq_t *seq = kseq_init(file);
    seq_counter sc = {.n_seq = 0, .sum_len = 0, .header_len=0};
    int l;
    // loop over sequences
    while ((l = kseq_read(seq)) >= 0) {
        if (sc.header_len == 0){
            sc.header_len = seq->name.l;
        }
        sc.n_seq++;
        sc.sum_len += (int) seq->seq.l;
    }
    return sc;
}



strdict_t *get_sequence_headers(const char * fname)
{
    // start by counting sequences and getting length of headers
    seq_counter sc;
    sc = count_sequences(fname);

    // fill set with sequence headers
    gzFile file = gzopen(fname, "r");
    kseq_t *seq = kseq_init(file);

    khash_t(str) *headers;
    headers = kh_init(str);
    khint_t k;
    kh_resize(str, headers, sc.n_seq);

    int l;
    int absent;
    int i = 0;
    char name[64];

    while ((l = kseq_read(seq)) >= 0) {
        strcpy(name, seq->name.s);
        k = kh_put(str, headers, name, &absent);
        assert(absent);
        kh_key(headers, k) = strdup(name);
        kh_val(headers, k) = i;
        i++;
    }
    gzclose(file);
    return headers;
}



void destroy_str_dict(strdict_t *h)
{
    khint_t k;
    for (k = 0; k < kh_end(h); ++k)
        if (kh_exist(h, k))
            free((char*)kh_key(h, k));
    kh_destroy(str, h);
}



strdict_t *get_index_headers(const char * fname)
{
    // open file stream
    FILE *fp = fopen(fname, "rb");
    // check idx validity and get n_seq
    char magic[4];
    uint32_t x[5];
    char name[64];
    if (fread(magic, 1, 4, fp) != 4) return 0;
    if (strncmp(magic, MM_IDX_MAGIC, 4) != 0) return 0;
    if (fread(x, 4, 5, fp) != 5) return 0;
    int n_seq = x[3];
    // init str set for headers
    khash_t(str) *idx_headers;
    idx_headers = kh_init(str);
    khint_t k = 0;
    kh_resize(str, idx_headers, n_seq);
    int absent;
    // iterate sequence in the idx
    for (int i = 0; i < n_seq; ++i) {
        uint8_t l;
        uint32_t seqlen;
        fread(&l, 1, 1, fp);
        if (l) {
            fread(name, 1, l, fp);
            k = kh_put(str, idx_headers, name, &absent);
            assert(absent);
            kh_key(idx_headers, k) = strdup(name);
            kh_val(idx_headers, k) = i;
        }
        // discard seq length
        fread(&seqlen, 4, 1, fp);
    }
    fclose(fp);
    return idx_headers;
}



intset_t *rel_comp(strdict_t *a, strdict_t *b)
{
    int_arr res;
    int r = 0;
    res.a = malloc(sizeof(khint_t) * kh_size(a));
    khint_t k;
    for (k = 0; k < kh_end(a); k++){
        // key does not exist in a
        if (!kh_exist(a, k)) continue;
        // key exists in a
        khint_t q;
        int is_missing;
        q = kh_get(str, b, kh_key(a, k));   // query b for the key at k in a
        is_missing = (q == kh_end(b));
        // key does not exist in b
        if (is_missing){
            res.a[r] = kh_val(a, k);
            r++;
        }
    }
    res.l = r;
    // put into int set
    khash_t(int_set) *h;
    h = kh_init(int_set);
    kh_resize(int_set, h, kh_size(a));  // needs to be same size for indexing later
    int absent;

    for (int i = 0; i < res.l; ++i) {
        k = kh_put(int_set, h, res.a[i], &absent);
        assert(absent);
    }
    return h;
}



Sketch sketch_sequence(kseq_t *seq, int w, int k, int r)
{
    Sketch sk = {.v.n = 0, .v.a = 0};
    strcpy(sk.n, seq->name.s);
    int l = (int) seq->seq.l;
    // definition is in mm2's sketch.c
    mm_sketch(0, seq->seq.s, l, w, k, r, 0, &sk.v);
    return sk;
}



Sketch *collect_sketches(const char * fname, intset_t *indices_add, int_arr *free_indices)
{
    // collect the sketches of sequences not in idx
    gzFile file = gzopen(fname, "r");
    kseq_t *seq = kseq_init(file);
    int nseq = kh_size(indices_add);

    // allocate mem for array of structs
    Sketch sk;
    Sketch *data = malloc(sizeof(Sketch) * nseq);

    int l;
    int i = 0;
    int j = 0;

    while ((l = kseq_read(seq)) >= 0) {
        if (!in_int_set(indices_add, i)) {
            i++;
            continue;
        }
        strcpy(data[j].n, seq->name.s);
        sk = sketch_sequence(seq, WINDOW_SIZE, KMER_SIZE, free_indices->a[j]);
        data[j].v.n = sk.v.n;
        data[j].v.m = sk.v.m;
        data[j].v.a = malloc(sizeof(mm128_t) * sk.v.n);
        data[j].v.a = sk.v.a;
        i++;
        j++;
    }
    gzclose(file);
    return data;
}



int_arr prep_indices_to_populate(intset_t *empty_indices, int n_seqs_in_idx, int n_new_seqs)
{
    int_arr res;
    int r = 0;

    for (int i = 0; i < n_seqs_in_idx; i++){
        if (!in_int_set(empty_indices, i)) continue;
        res.a[r] = kh_key(empty_indices, i);
        r++;
    }
    int k = 0;
    while (r < n_new_seqs){
        res.a[r] = n_seqs_in_idx + k;
        r++;
        k++;
    }
    res.l = r;
    return res;
}



void print_args(int argc, char *argv[]){
    // check correct number of args
    if (argc != 4) {
        printf("Usage: %s <in.fastx> <in.mmi> <out.mmi>\n", argv[0]);
        assert(0);
    }

    // print some info
    printf("mm2-ii %s %s \n", MM2_II_VERSION, MM2_VERSION);
    for (int i = 0; i < argc; i++) {
        printf("%s ", argv[i]);
    }
    printf("\n\n");
}



int main(int argc, char *argv[]) {
    // print usage and info
    print_args(argc, argv);

    // assign CL args
    const char * sequencefile = argv[1];
    const char * idx_in = argv[2];
    const char * idx_out = argv[3];

    // general timer
    clock_t tic = clock();

    // load headers of seq and idx files into sets
    clock_t tic01 = clock();
    khash_t(str) *seq_headers;
    khash_t(str) *idx_headers;
    seq_headers = get_sequence_headers(sequencefile);
    idx_headers = get_index_headers(idx_in);
    clock_t toc01 = clock();

    // find the difference between the index and the sequence file
    clock_t tic02 = clock();
    intset_t *a_not_b;
    intset_t *b_not_a;
    intset_t *idx_del;
    a_not_b = rel_comp(seq_headers, idx_headers);
    b_not_a = rel_comp(idx_headers, seq_headers);
    printf("seq to add: %d\n", kh_size(a_not_b));
    printf("seq to remove: %d\n", kh_size(b_not_a));
    idx_del = b_not_a;
    clock_t toc02 = clock();

    // prep indices to place new seqs   TODO
    //    int_arr free_indices;
    //    free_indices = prep_indices_to_populate(b_not_a, kh_size(idx_headers), kh_size(a_not_b));

    // collect new sketches, inject read ids
    //    Sketch *sk;
    //    sk = collect_sketches(sequencefile, a_not_b, &free_indices);


    // TODO add new sketches without having to resort buckets
    // loop sketches -> loop mms of sketch -> index into bucket ->

    // load and modify index
    clock_t tic03 = clock();
    FILE *fp_in = fopen(idx_in, "rb");
    FILE *fp_out = fopen(idx_out, "wb");

    // dump basics
    idx_info mm_info;
    dump_basics(fp_in, fp_out, &mm_info);

    // dump seq info
    dump_sequence_info(fp_in, fp_out, b_not_a, mm_info.x[3]);


    // process buckets
    for (int i = 0; i < 1<<mm_info.x[2]; ++i) {
        uint32_t size;
        int32_t n;   // size of the _p_ array
        uint64_t *p; // position array for minimizers appearing >1 times

        // read bucket basics
        fread(&n, 4, 1, fp_in);
        p = (uint64_t*)malloc(n * 8);
        fread(p, 8, n, fp_in);
        fread(&size, 4, 1, fp_in);

        // dump empty bucket
        if (size == 0){
            fwrite(&n, 4, 1, fp_out);
            fwrite(p, 8, n, fp_out);
            fwrite(&size, 4, 1, fp_out);
            continue;
        }

        // non-empty bucket
        idx_bucket b = {.s = size, .n = n, .p = p};
        b.keys = malloc(sizeof(uint64_t) * size);
        b.vals = malloc(sizeof(uint64_t) * size);
        uint64_t mm[2];

        // read minimizers in this bucket
        for (uint32_t j = 0; j < size; ++j) {
            fread(mm, 8, 2, fp_in);
            b.keys[j] = mm[0];
            b.vals[j] = mm[1];
        }


        // in-place modification of minimizers in the hash tables
        uint32_t n_del_s = 0;
        uint32_t n_del_m = 0;
        uint32_t n_red = 0;
        uint32_t n_mod = 0;
        int *pdel;


        // handle multitons and return indicator array for which elements to delete in p
        pdel = process_bucket(&b, idx_del, &n_del_s, &n_del_m, &n_red, &n_mod);
        // modify p array according to pdel array
        modify_parray(fp_out, &b, pdel);
        // write the modified minimmizers to file
        dump_minimizers(fp_out, &b, &n_del_s, &n_del_m);

        free(b.keys);
        free(b.vals);
        free(b.p);
        free(pdel);
    } // END OF BUCKET

    fflush(fp_out);
    fclose(fp_in);
    fclose(fp_out);

    clock_t toc03 = clock();

    printf("idx written to: %s\n", idx_out);

    // clean-up
    destroy_str_dict(seq_headers);
    destroy_str_dict(idx_headers);

    // timing analysis
    clock_t toc = clock();
    double tt = toc - tic;
    printf("t: %f s\n\n", tt / CLOCKS_PER_SEC);
    double t01 = (toc01 - tic01) / tt * 100;
    double t02 = (toc02 - tic02) / tt * 100;
    double t03 = (toc03 - tic03) / tt * 100;
    printf("01:\t %f \n", t01);
    printf("02:\t %f \n", t02);
    printf("03:\t %f \n", t03);

    return 0;
}


