#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>    // bool types
#include <sys/stat.h>   // stat_record, file checking

// external dependencies
#include <zlib.h>
#include <hashmap.h>

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

// magic identifier for mmis
#define MM_IDX_MAGIC    "MMI\2"
#define MM2_VERSION     "2.24-r1150"
#define MM2_II_VERSION  "0.0.1"
#define WINDOW_SIZE     5
#define KMER_SIZE       15
#define NTHREADS        1


typedef HASHMAP(char, mm128_v) sketchmap_t;


// main struct to hold minimizers & sequence name
typedef struct{
    char n[64];
    mm128_v v;
} Sketch;


typedef struct{
    int n, l;
} seq_counter;


// declarations for sketch-file IO
bool file_has_content(const char * fname);
bool write_sketches(const char *filename, Sketch *data, int total);
Sketch *read_sketches(const char *filename, int *total);
// put all sketches from file into a hashmap
bool fill_hashmap(sketchmap_t *hm, Sketch *data, const int *total);
// generate minimizer sketch from a single sequence
Sketch sketch_sequence(kseq_t *seq, int w, int k);
// merge arrays of loaded and new sketches for writing to file
Sketch *mergearray(Sketch *a, Sketch *b, int size_a, int size_b);
// count number of sequences in file and total length
seq_counter count_sequences(const char * fname);
// write mmi to file, wraps mm_idx_dump
bool write_index(const char * fname, mm_idx_t * mi);
// iterate over sequences and either put loaded or generated sketches into idx buckets
Sketch *fill_idx_buckets(const char * fname, mm_idx_t * mi, int n, sketchmap_t *hm, int *new_seqs);




int main(int argc, char *argv[])
{
    // check correct number of args
    if (argc != 4) {
        printf("Usage: %s <in.sketches> <in.fasta> <out.mmi>\n", argv[0]);
        return 1;
    }

    // print some info
    printf("mm2-ii %s %s \n", MM2_II_VERSION, MM2_VERSION);
    for (int i = 0; i < argc; i++) {
        printf("%s ", argv[i]);
    }
    printf("\n\n");

    const char * sketchfile = argv[1];
    const char * sequencefile = argv[2];
    const char * idxfile = argv[3];

    int total = 0;

    // init hashmap
    sketchmap_t hm;
    hashmap_init(&hm, hashmap_hash_string, strcmp);
    Sketch *sketches_read;

    if (file_has_content(sketchfile)){
        // read sketches from file
        sketches_read = read_sketches(sketchfile, &total);

        if (sketches_read == NULL){
            printf("Error reading sketchfile.\n");
            return 1;
        }

        printf("%d sketches loaded\n", total);

        // fill hashmap
        if (fill_hashmap(&hm, sketches_read, &total) != true){
            printf("Error filling hashmap\n");
            return 1;
        }
    } else{
        // dummy allocation
        sketches_read = malloc(sizeof(Sketch) * 1);
    }
    printf("hashmap created\n");

    // count input sequences and length
    seq_counter seq_count;
    seq_count = count_sequences(sequencefile);
    printf("%d sequences; %d total len\n", seq_count.n, seq_count.l);


    // initialise index - code mostly from mm_idx_str
    mm_idx_t * mi;
    mi = init_index(WINDOW_SIZE, KMER_SIZE);
    mi->n_seq = seq_count.n;
    //	mi->seq = (mm_idx_seq_t*)kcalloc(mi->km, nseq, sizeof(mm_idx_seq_t)); // ->seq is allocated from km
    mm_idx_seq_t * seqs = malloc(sizeof(mm_idx_seq_t) * seq_count.n);
    mi->seq = seqs;
    mi->S = (uint32_t*)calloc((seq_count.l + 7) / 8, 4);
    khash_t(str) *h;
    h = kh_init(str);
    mi->h = h;


    // iterate sequences and create index
    int new_seqs = 0;
    Sketch *sketches_new;
    sketches_new = fill_idx_buckets(sequencefile, mi, seq_count.n, &hm, &new_seqs);

    mm_idx_post(mi, NTHREADS);

    // save index to file
    if (write_index(idxfile, mi) != true){
        printf("Error writing index file\n");
    }
    printf("idx written to: %s\n", idxfile);


    // merge old and new sketches
    Sketch *merged_sketches;
    merged_sketches = mergearray(sketches_read, sketches_new, total, new_seqs);
    total += new_seqs;
    printf("%d total sketches; %d new sketches\n", total, new_seqs);

    // write new sketch file
    if (write_sketches(sketchfile, merged_sketches, total)){
        printf("%d Sketches written to: %s.\n", total, sketchfile);
    }
    else{
        printf("Error writing sketches.\n");
        return 1;
    }


    // clean-up
    hashmap_cleanup(&hm);
    free(seqs);
    free(sketches_read);
    free(sketches_new);
    free(merged_sketches);
    mm_idx_destroy(mi);
    return 0;
}




// return true iff file not empty
bool file_has_content(const char * fname) {
    struct stat stat_record;
    if (stat(fname, &stat_record)) {
        printf("File does not exist\n");
        return false;
    } else if (stat_record.st_size <= 1) {
        printf("File exists, is empty\n");
        return false;
    } else {
        printf("File exists, has content\n");
        return true;
    }
}



bool write_sketches(const char *filename, Sketch *data, int total)
{
    // total
    // data[0] name mm128_v.n mm128_v.m *mm128_t (minimizer values x, y)
    // data[1] ...
    FILE *file = fopen(filename, "wb");
    if (file == NULL) return false;
    // total number of structs
    if (fwrite(&total, sizeof(int), 1, file) != 1)
        return false;

    // iterate sketches, write mm arrays to file
    for (int t = 0; t < total; t++){
        fwrite(&data[t].n, sizeof(char[64]), 1, file);
        fwrite(&data[t].v.n, sizeof(size_t), 1, file);
        fwrite(&data[t].v.m, sizeof(size_t), 1, file);
        // write array values
        for (int u = 0; u < (int) data[t].v.n; u++){
            fwrite(&data[t].v.a[u], sizeof(mm128_t), 1, file);
        }
    }
    if (fclose(file) == EOF) return false;
    return true;
}



Sketch *read_sketches(const char *filename, int *total)
{
    // read structs from file, returns pointer to array
    // also returns number of structs by pointer to total
    FILE *file = fopen(filename, "rb");
    if (file == NULL) return NULL;

    // read total into total pointer
    if (fread(total, sizeof(int), 1, file) != 1)
        return NULL;

    // allocate mem for array of structs
    Sketch *data = malloc(sizeof(Sketch) * *total);

    for (int t = 0; t < *total; t++){
        char name[64];
        size_t n = 0, m = 0;
        fread(name, sizeof(char[64]), 1, file);
        fread(&n, sizeof(size_t), 1, file);
        fread(&m, sizeof(size_t), 1, file);
        strcpy(data[t].n, name);
        data[t].v.n = n;
        data[t].v.m = m;
        // read minimizer array for this sketch
        mm128_t * arr = malloc(sizeof(mm128_t) * n);
        fread(arr, sizeof(mm128_t), n, file);
        data[t].v.a = arr;
        free(arr);
    }
    fclose(file);
    return data; // pointer to arr of structs
}



bool fill_hashmap(sketchmap_t *hm, Sketch *data, const int *total){
    for (int m = 0; m < *total; m++) {
        char *n;
        mm128_v *v;
        n = strdup(data[m].n);
        v = &data[m].v;
        int insertion = hashmap_put(hm, n, v);
        if (insertion < 0) {
            printf("insertion failed: %s\n", strerror(-insertion));
            return false;
        }
    }
    return true;
}



Sketch sketch_sequence(kseq_t *seq, int w, int k)
{
    Sketch sk = {.v.n = 0, .v.a = 0};
    strcpy(sk.n, seq->name.s);
    int l = (int) seq->seq.l;
    // definition is in mm2's sketch.c
    mm_sketch(0, seq->seq.s, l, w, k, 0, 0, &sk.v);
    return sk;
}



Sketch *mergearray(Sketch *a, Sketch *b, int size_a, int size_b)
{
    int size_total = size_a + size_b; // calc final size
    int i, j;
    Sketch *c;
    c = malloc(sizeof(Sketch) * size_total);

    // copy sketches from array a into c
    for (i = 0; i < size_a; i++) {
        c[i] = a[i];
    }

    // copy sketches from array b into c
    for (i = 0, j = size_a;
         j < size_total && i < size_b; i++, j++) {
        c[j] = b[i];
    }
    return c;
}



seq_counter count_sequences(const char * fname){
    gzFile file = gzopen(fname, "r");
    kseq_t *seq = kseq_init(file);
    seq_counter sc = {.n = 0, .l = 0};
    int l;

    while ((l = kseq_read(seq)) >= 0) {
        sc.n++;
        sc.l += (int) seq->seq.l;
    }
    kseq_destroy(seq);
    gzclose(file);
    return sc;
}



bool write_index(const char * fname, mm_idx_t * mi){
    FILE *file = fopen(fname, "wb");
    if (file == NULL) return false;
    mm_idx_dump(file, mi);
    if (fclose(file) == EOF) return false;
    return true;
}



Sketch *fill_idx_buckets(const char * fname, mm_idx_t * mi, int n, sketchmap_t *hm, int *new_seqs){
    // fill index buckets and return array of new sketches
    gzFile file = gzopen(fname, "r");
    kseq_t *seq = kseq_init(file);
    int l;
    int idx = 0;
    int offset_c = 0;
    mm128_v *q;                 // hashmap value
    Sketch sketch_new;          // sketch of new sequence
    Sketch *sketches_new;       // array of new sketches
    sketches_new = malloc(sizeof(Sketch) * n);

    // code mostly equivalent to mm_idx_str
    while ((l = kseq_read(seq)) >= 0) {
        const char *s = seq->seq.s;
        mm_idx_seq_t *p = &mi->seq[idx];
        uint32_t j;
        int absent;
        // p->name = (char*)kmalloc(mi->km, strlen(name[i]) + 1);
        p->name = malloc(strlen(seq->name.s) + 1);
        strcpy(p->name, seq->name.s);
        kh_put(str, mi->h, p->name, &absent);
        assert(absent);

        p->offset = offset_c;
        p->len = seq->seq.l;
        p->is_alt = 0;
        for (j = 0; j < p->len; ++j) {
            int c = seq_nt4_table[(uint8_t)s[j]];
            uint64_t o = offset_c + j;
            mm_seq4_set(mi->S, o, c);
        }
        int plen; plen = (int) p->len;
        offset_c += plen;

        // either find sketch in hashmap or create
        q = hashmap_get(hm, seq->name.s);
        if (!q){
            // printf("%s NO\n", seq->name.s);
            sketch_new = sketch_sequence(seq, WINDOW_SIZE, KMER_SIZE);
            sketches_new[*new_seqs] = sketch_new;
            (*new_seqs)++;
            int sn = (int) sketch_new.v.n;
            mm_idx_add(mi, sn, sketch_new.v.a);
        }
        else {
            // printf("%s YES\n", seq->name.s);
            int qn = (int) q->n;
            mm_idx_add(mi, qn, q->a);
        }
        idx++;
    }
    kseq_destroy(seq);
    gzclose(file);
    return sketches_new;
}



