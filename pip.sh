#!/bin/bash
script='/u/home/y/yangao/program/NanoClu/NanoClu.py'
file=(/u/home/y/yangao/data/snyder_pnas/GM12878.fa.gz ~/data/nanopore/Albacore_20170317/2d.fq ~/data/nanopore/Albacore_20170427/2d.fq ~/data/nanopore/combine/com_2d.fq)
file_N=4

output=(/u/home/y/yangao/data/pb/GM12878/new_base_on_GMAP/shared_splice_site /u/home/y/yangao/data/nanopore/Albacore_20170317/new_base_on_GMAP/shared_splice_site /u/home/y/yangao/data/nanopore/Albacore_20170427/new_base_on_GMAP/shared_splice_site ~/data/nanopore/combine/shared_splice_site)

mode=(0)

bin_min_cnt=(2)

clu_min_size=(5)

thread_n=16

for ((i=0; i<$file_N; ++i))
do
    fq=${file[$i]}
    out=${output[$i]}
    for m in "${mode[@]}"
    do
        for b in "${bin_min_cnt[@]}"
        do
            for c in "${clu_min_size[@]}"
            do
                echo "python $script -f $fq -d hg19_M -t $thread_n -o $out -M $m --min-bin-cnt $b --min-clu-size $c -u UPDATE_GTF  -p ~/software/poaV2/poa"
                #python $script -f $fq -d hg19_M -t $thread_n -o $out -M $m --min-bin-cnt $b --min-clu-size $c  -u UPDATE_GTF  -p ~/software/poaV2/poa
            done
        done
    done
done
