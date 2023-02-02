# mm2-ii
Incremental indices for minimap2


## General idea

This tool speeds up alignment of reads by avoiding to re-index sequences that have been indexed before.
E.g. if all-vs-all alignments of reads are performed and new reads shall be added to the comparison,
then indexing reads that were already indexed previously would be wasteful.


## More detailed procedure

TODO


## TODOs

- currently only removal of sequences works
- adding new sequences is work in progress


## Installation & Usage


clone this repository and also clone [minimap2](https://github.com/lh3/minimap2) (currently tested with 2.24-r1150) into it (needed for compiler):

```
git clone https://github.com/W-L/mm2-ii.git && cd mm2-ii && git clone https://github.com/lh3/minimap2.git
```

create build directory and run the build:

```
mkdir -p build && cmake -S . -B build && cd build && make
```

The executable should be available as: `./build/mm2ii`

## Usage

`mm2ii <in.fastx> <in.mmi> <out.mmi>`

- `<in.fastx>`: sequence file for which an index should be created.
- `<in.mmi>`: input index to be modified.
- `<out.mmi>`: output name of modified index.


## License

Licensed under MIT






