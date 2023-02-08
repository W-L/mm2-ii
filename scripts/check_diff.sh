#DATA=$1
#DATA0="data/zymolive_00010.short.fa"
#DATA1="data/zymolive_00010m3.short.fa"

#DATA0="data/zymolive_01000.fa"
#DATA1="data/zymolive_00100.fa"

#DATA0="data/zymo_50k.fq"
#DATA1="data/zymo_49999.fq"

DATA0=$1
DATA1=$2



D0=out/idx_0.mmi
D1=out/idx_1.mmi
D2=out/idx_1_custom.mmi
R0=out/z_orig.paf
R1=out/z_custom.paf

rm $D0 $D1 $D2 $R0 $R1


## create two indices: one with full set, one with reduced set
echo ""
minimap2 -x map-ont -k15 -w5 -d $D0 "$DATA0"
echo ""
minimap2 -x map-ont -k15 -w5 -d $D1 "$DATA1"



## run mm2ii
printf '\n\n==========\n\n'
./src/cmake-build-debug/mm2ii $DATA1 $D0 $D2 || rm $D2   ## remove partial index
#./src/cmake-build-release/mm2ii $DATA1 $D0 $D2 || rm $D2   ## remove partial index
printf '\n\n==========\n\n'
## check diff between original and custom
#diff $D1 $D2 && printf "\n --- no diff --- \n \n"



## map with original and custom indices
echo "reduced index on full data"
minimap2 -x map-ont -k15 -w5 $D1 "$DATA0" >$R0
echo "custom index on full data"
minimap2 -x map-ont -k15 -w5 $D2 "$DATA0" >$R1
echo ""

## check if paf files are different
diff <(cut -f1-9 $R0)  <(cut -f1-9 $R1) >idxdiff
ls -lh idxdiff
#diff $R0 $R1 && printf "\n --- no diff --- \n \n"

## additional checks for binary and hex representation of indices
## not very useful if sequences retain their IDs compared to newly created ones
#xxd -b $D1 >$D1.xxd
#xxd -b $D2 >$D2.xxd
#diff $D1.xxd $D2.xxd >idxdiff
