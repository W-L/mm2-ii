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


## Dependencies, installation & Usage

Dependencies are:

- [hashmap](https://github.com/DavidLeeds/hashmap)

```
git clone https://github.com/DavidLeeds/hashmap.git
mkdir build-hashmap && cd build-hashmap
cmake ../hashmap && make
sudo make install
```


- [minimap2](https://github.com/lh3/minimap2) (currently tested with 2.24-r1150)

(does not need to be installed, just present at compile-time)


Installation:

clone this repository and also clone minimap2 into it (needed for compiler):

`git clone https://github.com/W-L/mm2-ii.git && cd mm2-ii && git clone https://github.com/lh3/minimap2.git`

create build directory and run the build:

`mkdir -p build && cmake -S . -B build && cd build && make`

The executable should then be available as: `./build/mm2ii`

## Usage

`mm2ii <in.sketches> <out.sketches> <in.fastx> <out.mmi>`

- `<in.sketches>`: file containing minimizer sketches, can also not exist or be empty.
- `<out.sketches>`: filename for writing minimizer sketched to.
- `<in.fastx>`: sequence file for which an index should be created. 
Minimizers only need to be calculated for sequences that do not already have a sketch in `sketch.in`.
- `<out.mmi>`: output name of the minimap2 index.




## License

Licensed under MIT






