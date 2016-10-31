export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/youyir/src/anaconda/pkgs/gcc-4.8.3-1/lib
../measure_adj/measure_adj
#./prepare_adj_src.pl -m CMTSOLUTION  -s STATIONS -i OUTPUT_FILES/*.adj
mv OUTPUT_FILES/*.adj OUTPUT_FILES/KBL.IU.60_100/
