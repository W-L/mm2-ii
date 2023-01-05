# mm2-ii
Incremental indices for minimap2


## General idea

This tool speeds up alignment of reads by avoiding to re-index sequences that have been indexed before.
E.g. if all-vs-all alignments of reads are performed and new reads shall be added to the comparison,
then indexing reads that were already indexed previously would be wasteful.

Here, we create an intermediate file that holds the minimizer sketches of a collection of sequences (fasta/q).
We load a file that contains such sketches (if present) and a sequence file. 
A minimap2 index is then created for all sequences contained in the sequence file. 
However, we only need to find minimizer sketches for sequences that are not contained in the sketch-file.
We then generate the minimap2 index and write all minimizer sketches to a file

## More detailed procedure

- Initialise a hashmap for minimizer sketches
- Load sketches from file and fill the hashmap; or continue with empty hashmap
- Get number of sequences and total sequence length from sequence file (needed to initialise idx)
- Initialise idx (incl. khash, seqs, ...)
- Fill idx buckets with minimizer sketches from sequences
- Run idx post and write idx to file
- Merge loaded and created minimizer sketches and write them to sketch file
- Perform clean-up


## Installation & dependencies

Dependencies are:

- [minimap2](https://github.com/lh3/minimap2)
- [hashmap](https://github.com/DavidLeeds/hashmap)


Installation:

`git clone TODO`

`mkdir -p build && cmake -S . -B build && cd build && make`

Executable should be available as: `./build/mm2ii`

## Usage

`mm2ii <in.sketches> <in.fasta> <out.mmi>`


## TODO

- the index is not completely identical to an index produced by mm2



