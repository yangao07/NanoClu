#!/usr/bin/env bash
poa=$1
fxtools=$2
clu_in=$3
clu_out_tmp=$4
poa_mat=$5
base_name=$6
clu_out=$7
#clu_in_n=$8

# TODO: -do_progressive -best
$poa -read_fasta $clu_in -pir $clu_out_tmp  -tolower -best -hb $poa_mat 2> /dev/null
sed -i "s/>CON/>${base_name}_CON/" $clu_out_tmp
$fxtools fn -m consensus $clu_out_tmp > $clu_out
rm $clu_out_tmp