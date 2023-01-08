import mappy
import sys

idx = sys.argv[1]
fa = sys.argv[2]


m = mappy.Aligner(fn_idx_in=idx, preset="map-ont", k=15, w=5)

for n in m.seq_names:
    print(n)

print()

for name, seq, qual in mappy.fastx_read(fa):
    for hit in m.map(seq):
        print("{}\t{}\t{}".format(hit.ctg, hit.r_st, hit.r_en))

print()
