#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>    // bool types
#include <sys/stat.h>   // stat_record, file checking
#include <inttypes.h>

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
#define MM2_VERSION     "2.24-r1150"
#define MM2_II_VERSION  "0.0.1"
#define WINDOW_SIZE     5
#define KMER_SIZE       15
#define NTHREADS        1


typedef HASHMAP(char, mm128_t) vmap_t;
typedef HASHMAP(char, size_t) sizemap_t;


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
bool fill_hashmap(vmap_t *vm, sizemap_t *sm, Sketch *data, const int *loaded);
// generate minimizer sketch from a single sequence
Sketch sketch_sequence(kseq_t *seq, int w, int k, int r);
// count number of sequences in file and total length
seq_counter count_sequences(const char * fname);
// write mmi to file, wraps mm_idx_dump
bool write_index(const char * fname, mm_idx_t * mi);
// collect the loaded and generated sketches in an array
void collect_sketches(const char * fname, mm_idx_t * mi, vmap_t *vm, sizemap_t *sm, Sketch * cs);
// iterate over sequences and put into idx buckets
void fill_idx_buckets(mm_idx_t * mi, int n, Sketch *data);
// modifying the y-values in mm sketches is necessary
uint64_t modify_y(uint64_t y, int rid);
void change_rid(mm128_t *arr, size_t n, int rid);
// helper func for debugging
void print_sketches(Sketch *data, int n);





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


    // count input sequences and length
    seq_counter seq_count;
    seq_count = count_sequences(sequencefile);
    printf("%d sequences; %d total len\n", seq_count.n, seq_count.l);


    // init hashmap
    vmap_t vm;
    sizemap_t sm;
    hashmap_init(&vm, hashmap_hash_string, strcmp);
    hashmap_init(&sm, hashmap_hash_string, strcmp);
    Sketch *sketches_read;

    int loaded = 0;

    if (file_has_content(sketchfile)){
        // read sketches from file
        sketches_read = read_sketches(sketchfile, &loaded);

        // print_sketches(sketches_read, loaded);

        if (sketches_read == NULL){
            printf("Error reading sketchfile.\n");
            return 1;
        }

        printf("%d sketches loaded\n", loaded);

    } else{
        // dummy allocation
        sketches_read = malloc(sizeof(Sketch) * 1);
    }


    // initialise index - code mostly from mm_idx_str
    mm_idx_t * mi;
    mi = init_index(WINDOW_SIZE, KMER_SIZE);
    mi->n_seq = seq_count.n;
    mm_idx_seq_t * seqs = malloc(sizeof(mm_idx_seq_t) * seq_count.n);
    mi->seq = seqs;
    mi->S = (uint32_t*)calloc((seq_count.l + 7) / 8, 4);
    khash_t(str) *h;
    h = kh_init(str);
    mi->h = h;


    // fill hashmap
    if (fill_hashmap(&vm, &sm, sketches_read, &loaded) != true){
        printf("Error filling hashmap\n");
        return 1;
    }
    printf("hashmap created\n");

    // iterate sequences and create index
    printf("collecting sketches\n");
    Sketch *collected_sketches;
    collected_sketches = malloc(seq_count.n * sizeof(mm128_t) * 20000); // TODO
    collect_sketches(sequencefile, mi, &vm, &sm, collected_sketches);

    // print_sketches(collected_sketches, seq_count.n);

    // add sequences to index and run post-processing
    fill_idx_buckets(mi, seq_count.n, collected_sketches);
    mm_idx_post(mi, NTHREADS);

    // save index to file
    if (write_index(idxfile, mi) != true){
        printf("Error writing index file\n");
        return 1;
    }
    printf("idx written to: %s\n", idxfile);


    // write new sketch file
    if (write_sketches(sketchfile, collected_sketches, seq_count.n)){
        printf("%d Sketches written to: %s.\n", seq_count.n, sketchfile);
    }
    else{
        printf("Error writing sketches.\n");
        return 1;
    }

    // clean-up
    hashmap_cleanup(&vm);
    hashmap_cleanup(&sm);
    free(seqs);
    free(sketches_read);
    free(collected_sketches);
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

        int array_size = (int) data[t].v.n;
        mm128_t * arr = malloc(sizeof(mm128_t) * array_size);

        for (int u = 0; u < array_size; u++){
            arr[u].x = data[t].v.a[u].x;
            arr[u].y = data[t].v.a[u].y;
        }

        // write the array of a single sketch to file
        fwrite(arr, sizeof(arr[0]), array_size, file);
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
        mm128_t *arr = malloc(sizeof(mm128_t) * n);
        fread(arr, sizeof(mm128_t), n, file);
        data[t].v.a = arr;
    }
    fclose(file);
    return data; // pointer to arr of structs
}



bool fill_hashmap(vmap_t *hm, sizemap_t *sm, Sketch *data, const int *loaded){
    for (int m = 0; m < *loaded; m++) {
        size_t *array_size = &data[m].v.n;
        mm128_t *arr = malloc(sizeof(mm128_t) * *array_size);

        for (int u = 0; u < *array_size; u++){
            arr[u].x = data[m].v.a[u].x;
            arr[u].y = data[m].v.a[u].y;
        }

        int ix = hashmap_put(hm, data[m].n, arr);
        if (ix < 0) {
            printf("insertion failed: %s\n", strerror(-ix));
            return false;
        }

        int is = hashmap_put(sm, data[m].n, array_size);
        if (is < 0) {
            printf("insertion failed: %s\n", strerror(-is));
            return false;
        }
    }
    return true;
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



seq_counter count_sequences(const char * fname){
    gzFile file = gzopen(fname, "r");
    kseq_t *seq = kseq_init(file);
    seq_counter sc = {.n = 0, .l = 0};
    int l;

    while ((l = kseq_read(seq)) >= 0) {
        sc.n++;
        sc.l += (int) seq->seq.l;
    }
    return sc;
}



bool write_index(const char * fname, mm_idx_t * mi){
    FILE *file = fopen(fname, "wb");
    if (file == NULL) return false;
    mm_idx_dump(file, mi);
    if (fclose(file) == EOF) return false;
    return true;
}



void fill_idx_buckets(mm_idx_t * mi, int n, Sketch *data){
    // fill index buckets and return array of new sketches
    printf("filling buckets\n");

    // code mostly equivalent to mm_idx_str
    for (int i = 0; i < n; i++){
        int sn = (int) data[i].v.n;
        mm_idx_add(mi, sn, data[i].v.a);
    }
}


uint64_t modify_y(uint64_t y, int rid){
    // take the y-value of a mm128_t as input and a rid to be injected
    // y is a combination of rid, position and strand
    // mm_sketch() in sketch.c: y = (uint64_t)rid<<32 | (uint32_t)i<<1 | z;
    // bitshifts and ORs
    uint64_t newy;
    uint32_t y_lowbits;
    // cast input to 32bit to lose the original rid
    y_lowbits = (uint32_t) y;
    newy = (uint64_t) rid<<32 | y_lowbits;
    // printf("%" PRIu64 "\n", y);
    return newy;
}


void change_rid(mm128_t *arr, size_t n, int rid){
    // iterate array and replace the rid encoded in y values
    for (int i = 0; i < n; i++){
        uint64_t newy;
        mm128_t *a = &arr[i];
        newy = modify_y(a->y, rid);
        a->y = newy;
    }
}




void collect_sketches(const char * fname, mm_idx_t * mi, vmap_t *vm, sizemap_t *sm, Sketch * cs){
    // fill index buckets and return array of new sketches
    gzFile file = gzopen(fname, "r");
    kseq_t *seq = kseq_init(file);
    int l;
    int idx = 0;
    int offset_c = 0;
    mm128_t *qv;                 // hashmap value
    int qs;

    // code mostly equivalent to mm_idx_str
    while ((l = kseq_read(seq)) >= 0) {
        const char *s = seq->seq.s;
        mm_idx_seq_t *p = &mi->seq[idx];
        uint32_t j;
        int absent;
        Sketch sketch_new;          // sketch of new sequence
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
        qv = hashmap_get(vm, seq->name.s);

        if (!qv){
            // printf("%s NO\n", seq->name.s);
            sketch_new = sketch_sequence(seq, WINDOW_SIZE, KMER_SIZE, idx);
            cs[idx] = sketch_new;
        }
        else {
            // printf("%s YES\n", seq->name.s);
            qs = (int) *hashmap_get(sm, seq->name.s);
            strcpy(sketch_new.n, seq->name.s);
            sketch_new.v.n = qs;
            sketch_new.v.m = 0;
            sketch_new.v.a = malloc(sizeof(mm128_t) * qs);

            for (int i = 0; i < qs; i++){
                sketch_new.v.a[i].x = qv[i].x;
                sketch_new.v.a[i].y = qv[i].y;
            }
            // inject the correct read id into loaded sketches
            change_rid(sketch_new.v.a, (size_t) qs, idx);

            cs[idx] = sketch_new;
        }
        idx++;
    }
}




void print_sketches(Sketch *data, int n){
    for (int i = 0; i < n; i++){
        int s = data[i].v.n;
        printf("%d\n", s);

        for (int j = 0; j < s; j++) {
            printf("%" PRIu64 " ", data[i].v.a[j].x);
        }
        printf("\n");
        for (int j = 0; j < s; j++) {
            printf("%" PRIu64 " ", data[i].v.a[j].y);
        }
        printf("\n\n\n");
    }
}
