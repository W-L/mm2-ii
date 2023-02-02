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

#include "index_mm2ii.h"

// inits for kseq & khash
KSEQ_INIT(gzFile, gzread)

#define idx_hash(a) ((a)>>1)
#define idx_eq(a, b) ((a)>>1 == (b)>>1)

KHASH_INIT(idx, uint64_t, uint64_t, 1, idx_hash, idx_eq)
typedef khash_t(idx) idxhash_t;
KHASH_MAP_INIT_STR(str, uint32_t)


KHASH_SET_INIT_INT(int_set)
typedef khash_t(int_set) intset_t;


KHASH_MAP_INIT_INT(int_dict, int)
typedef khash_t(int_dict) intdict_t;


//KHASH_MAP_INIT_STR(str_dict, int)
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



int in_int_set(intset_t *ih, const uint32_t query)
{
    khint_t q;
    int in_set;
    q = kh_get(int_set, ih, query);
    in_set = (q != kh_end(ih));
    if (in_set) return 1;
    return 0;
}


uint64_t dump_sequence_info(FILE *fp, const mm_idx_t *mi, intset_t *idx_del)
{
    uint64_t sum_len = 0;
    uint32_t i;
    for (i = 0; i < mi->n_seq; ++i) {
        if (mi->seq[i].name) {
            uint8_t l = strlen(mi->seq[i].name);
            fwrite(&l, 1, 1, fp);
            fwrite(mi->seq[i].name, 1, l, fp);
        } else {
            uint8_t l = 0;
            fwrite(&l, 1, 1, fp);
        }

        // for deleted seqs: set seq length to 0
        // not deleting sequences, to avoid shuffling y-values
        if (in_int_set(idx_del, i)){
            mi->seq[i].len = 0;
        }
        fwrite(&mi->seq[i].len, 4, 1, fp);
        sum_len += mi->seq[i].len;
    }
    return sum_len;
}



void dump_basics(FILE *fp, const mm_idx_t *mi)
{
    uint32_t x[5];
    x[0] = mi->w, x[1] = mi->k, x[2] = mi->b, x[3] = mi->n_seq, x[4] = mi->flag;
    fwrite(MM_IDX_MAGIC, 1, 4, fp);
    fwrite(x, 4, 5, fp);
}





void dump_without_mod(FILE *fp, mm_idx_bucket_t *b)
{
    idxhash_t *h = (idxhash_t*)b->h;
    uint32_t size = h? h->size : 0;
    khint_t k;

    fwrite(&size, 4, 1, fp);
    for (k = 0; k < kh_end(h); ++k) {
        uint64_t mm[2];
        if (!kh_exist(h, k)) continue;
        mm[0] = kh_key(h, k), mm[1] = kh_val(h, k);
        fwrite(mm, 8, 2, fp);
    }
}




void reduce_to_singleton(mm_idx_bucket_t *b, int ppos, uint64_t mnm)
{
    idxhash_t  *h = (idxhash_t *)b->h;
    khint_t k_ind;
    k_ind = kh_get(idx, h, mnm);
    uint64_t mm[2];
    mm[0] = kh_key(h, k_ind), mm[1] = kh_val(h, k_ind);
    // new minimizer
    uint64_t nn[2];
    nn[0] = mm[0]|1, nn[1] = b->p[ppos];
    // delete old key
    kh_del(idx, h, k_ind);
    // assign new key value in hash table
    int absent;
    khint_t k;
    k = kh_put(idx, h, nn[0], &absent);
    assert(absent);
    kh_value(h, k) = nn[1];
}


void change_modified_multiton(mm_idx_bucket_t *b, int pshift, int newsize, uint64_t mnm){
    idxhash_t  *h = (idxhash_t *)b->h;
    khint_t k_ind;
    k_ind = kh_get(idx, h, mnm);
    uint64_t mm[2];
    mm[0] = kh_key(h, k_ind), mm[1] = kh_val(h, k_ind);
    uint32_t ppos = mm[1] >> 32;
    ppos = ppos - pshift;
    // modify the value of this mm
    kh_value(h, k_ind) = (uint64_t) ppos<<32 | newsize;
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


void modify_parray(FILE *fp, mm_idx_bucket_t *b, int *pdel)
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





int *process_multitons(mm_idx_bucket_t *b, intset_t *idx_del, int *n_del_multiton, int *n_del_reductions, int *n_mod){
    idxhash_t *h = (idxhash_t*)b->h;
    int m = b->n;
    int absent;
    khint_t k;

    // create a k-array with same length as p but pointing back to the hash table
    int mm_arr[m]; // = malloc(sizeof(int) * b->n);
    int cn_arr[m];
    int newsize_arr[m];
    int *pdel = malloc(sizeof(int) * m);

    // no multitons in this bucket
    if (m == 0){
        return pdel;
    }

    // collect multiton info arrays
    for (k = 0; k < kh_end(h); ++k) {
        if (!kh_exist(h, k)) continue;
        uint64_t mm[2];
        mm[0] = kh_key(h, k), mm[1] = kh_val(h, k);


        if ((mm[0] & 1) != 1) {   // multiton
            // get copy number
            uint32_t copyn;
            copyn = (uint32_t) mm[1];
            uint32_t ppos = mm[1] >> 32;
            int n_del_units = 0;

            // fill the respective positions in karr and carr
            for (int l = 0; l < copyn; l++) {
                mm_arr[ppos + l] = mm[0];
                cn_arr[ppos + l] = copyn;

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
        }
    }
    // effect multiton info arrays
    for (int l = 0; l < m;) {
        uint64_t mnm = mm_arr[l];
        int cn = cn_arr[l];

        // get new size
        int newsize = newsize_arr[l];
        assert(newsize >= 0);

        // no change
        if (cn == newsize){
            l += cn;
            continue;}

        // wipe-out
        if (newsize == 0){
            khint_t k_ind = kh_get(idx, h, mnm);
            kh_del(idx, h, k_ind);
            (*n_del_multiton)++;
            l += cn;
            continue;
        }

        // reduce to singleton
        if (newsize == 1){
            // search for the copy in p that comes from a non-deleted sequence
            int done = 0;
            for (int i = 0; i < cn; i++) {
                if (pdel[l + i] == 0) {
                    assert(!done);
                    reduce_to_singleton(b, l + i, mnm);
                    pdel[l + i] = 1;
                    done = 1;
                }
            }
            assert(done);
            (*n_del_reductions)++;
            l += cn;
            continue;
        }

        // remain multiton
        if (newsize > 1){
            // get pshift
            int *pdel_cumsum;
            pdel_cumsum = array_cumsum(pdel, m);
            int pshift = pdel_cumsum[l];
            change_modified_multiton(b, pshift, newsize, mnm);
            (*n_mod)++;
            l += cn;
            continue;
        }
    }
    return pdel;
}



void process_singletons(mm_idx_bucket_t *b, intset_t *idx_del, int *n_del_singleton)
{
    idxhash_t *h = (idxhash_t*)b->h;
    uint32_t size = h? h->size : 0;
    khint_t k;
    if (size == 0){return;}

    for (k = 0; k < kh_end(h); ++k) {
        if (!kh_exist(h, k)) continue;
        uint64_t mm[2];
        mm[0] = kh_key(h, k), mm[1] = kh_val(h, k);

        if ((mm[0] & 1) == 1) {   // singleton
            // check if source sequence is deleted
            if (in_int_set(idx_del, mm[1] >> 32)) {
                // delete mm
                (*n_del_singleton)++;
                kh_del(idx, h, k);
            }
        }
    }
}


void dump_hashtable(FILE *fp, mm_idx_bucket_t *b)
{
    idxhash_t *h = (idxhash_t*)b->h;
    uint32_t size = h? h->size : 0;
    khint_t k;

    // empty bucket
    if (size == 0){
        fwrite(&size, 4, 1, fp);
        return;
    }

    // non-empty bucket
    // size has alrady changed from dynamic updates to hashtable
    // no need to recalc
    fwrite(&size, 4, 1, fp);

    for (k = 0; k < kh_end(h); ++k) {
        if (!kh_exist(h, k)) continue;
        uint64_t mm[2];
        mm[0] = kh_key(h, k), mm[1] = kh_val(h, k);
        fwrite(mm, 8, 2, fp);
    }

}



void mm_idx_dump_mod(FILE *fp, const mm_idx_t *mi, intset_t *idx_del)
{
    // dump the magic and 5 basic parameters
    dump_basics(fp, mi);

    // dump sequence info
    uint64_t sum_len;
    sum_len = dump_sequence_info(fp, mi, idx_del);

    // iterate buckets for dumping
    uint32_t i;
    for (i = 0; i < 1<<mi->b; ++i) {
        mm_idx_bucket_t *b = &mi->B[i];
        idxhash_t  *h = (idxhash_t *)b->h;
        uint32_t size = h? h->size : 0;


        // this bucket does not have any mms in it
        // i.e. it does not have an initialised hash table
        if (size == 0){
            fwrite(&b->n, 4, 1, fp);
            fwrite(b->p, 8, b->n, fp);
            fwrite(&size, 4, 1, fp);
            continue;
        }


        // new approach: excplicit handling of multitons and singletons
        // in-place modification of minimizers in the hash tables
        int n_del_singleton = 0;
        int n_del_multiton = 0;
        int n_del_reductions = 0;
        int n_mod = 0;
        int *pdel;

        // handle multitons and return indicator array for which elements to delete in p
        pdel = process_multitons(b, idx_del, &n_del_multiton, &n_del_reductions, &n_mod);
        // delete singletons in-place
        process_singletons(b, idx_del, &n_del_singleton);
        // modify p array according to pdel array
        modify_parray(fp, b, pdel);
        // write the modified hash tables to file
        dump_hashtable(fp, b);
    }
    if (!(mi->flag & MM_I_NO_SEQ))
        fwrite(mi->S, 4, (sum_len + 7) / 8, fp);
    fflush(fp);
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

    while ((l = kseq_read(seq)) >= 0) {
        k = kh_put(str, headers, seq->name.s, &absent);
        assert(absent);
        kh_key(headers, k) = strdup(seq->name.s);
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

//
//mm_idx_t *mm_idx_init(int w, int k, int b, int flag)
//{
//    mm_idx_t *mi;
//    if (k*2 < b) b = k * 2;
//    if (w < 1) w = 1;
//    mi = (mm_idx_t*)calloc(1, sizeof(mm_idx_t));
//    mi->w = w, mi->k = k, mi->b = b, mi->flag = flag;
//    mi->B = (mm_idx_bucket_t*)calloc(1<<b, sizeof(mm_idx_bucket_t));
//    if (!(mm_dbg_flag & 1)) mi->km = km_init();
//    return mi;
//}



strdict_t *get_index_headers(const char * fname)
{
    // open file stream
    FILE *fp = fopen(fname, "rb");
    // check idx validity and get n_seq
    char magic[4];
    uint32_t x[5];
    mm_idx_t *mi;
    if (fread(magic, 1, 4, fp) != 4) return 0;
    if (strncmp(magic, MM_IDX_MAGIC, 4) != 0) return 0;
    if (fread(x, 4, 5, fp) != 5) return 0;
    mi = mm_idx_init(x[0], x[1], x[2], x[4]);

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
        char *name;
        fread(&l, 1, 1, fp);
        if (l) {
            name = (char*)kmalloc(mi->km, l + 1);
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


mm_idx_t *load_index(const char *idx_in){
    FILE *idx_read_fp = fopen(idx_in, "rb");
    mm_idx_t *mi;
    mi = mm_idx_load(idx_read_fp);
    fclose(idx_read_fp);
    return mi;
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
    khash_t(str) *seq_headers;
    khash_t(str) *idx_headers;
    seq_headers = get_sequence_headers(sequencefile);
    idx_headers = get_index_headers(idx_in);

    // find the difference between the index and the sequence file
    intset_t *a_not_b;
    intset_t *b_not_a;
    a_not_b = rel_comp(seq_headers, idx_headers);
    b_not_a = rel_comp(idx_headers, seq_headers);
    printf("seq to add: %d\n", kh_size(a_not_b));
    printf("seq to remove: %d\n", kh_size(b_not_a));

    // prep indices to place new seqs
    //    int_arr free_indices;
    //    free_indices = prep_indices_to_populate(b_not_a, kh_size(idx_headers), kh_size(a_not_b));

    // collect new sketches, inject read ids
    //    Sketch *sk;
    //    sk = collect_sketches(sequencefile, a_not_b, &free_indices);


    // TODO add new sketches without having to resort buckets
    // loop sketches -> loop mms of sketch -> index into bucket ->

    // load index from file
    mm_idx_t *mi = load_index(idx_in);


    // TODO move mod loop instead of during dumping


    // write index to file incl modifications
    FILE *idx_dump_fp = fopen(idx_out, "wb");
    mm_idx_dump_mod(idx_dump_fp, mi, b_not_a);
    fclose(idx_dump_fp);

    printf("idx written to: %s\n", idx_out);

    // clean-up
    mm_idx_destroy(mi);
    destroy_str_dict(seq_headers);
    destroy_str_dict(idx_headers);
    //    kh_destroy(int_set, h);


    // timing analysis
    clock_t toc = clock();
    double tt = toc - tic;
    printf("t: %f s\n\n", tt / CLOCKS_PER_SEC);
    //    double t0 = (toc0 - tic0) / tt * 100;
    //    printf("counting:\t %f \n", t001);

    return 0;
}


