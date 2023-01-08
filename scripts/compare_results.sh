
DATA=$1

./cmake-build-debug/mm2ii sketches.bin "$DATA" idx.mmi

python scripts/print_mmi_names.py idx.mmi "$DATA"

minimap2 -x map-ont -k 15 -w 5 -d idx_orig.mmi "$DATA" "$DATA" >z_orig.paf

minimap2 -x map-ont idx.mmi "$DATA" >z_custom.paf

echo
diff idx.mmi idx_orig.mmi
echo

diff z_orig.paf z_custom.paf

xxd -b idx.mmi >idx.mmi.xxd

xxd -b idx_orig.mmi >idx_orig.mmi.xxd

#diff idx.mmi.xxd idx_orig.mmi.xxd




