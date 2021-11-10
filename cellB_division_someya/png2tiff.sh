#!/bin/sh
#seqは連番を出力、%03gは全体３桁になるように０パディング（001-400）
for i in `seq -f %03g 1 691`

#png→tiff
do
    echo $i
    for j in `seq -f %01g 0 12`
    do
        convert "mark${j}/p${i}.png" "mark${j}/t${i}.tiff"
    done
done

#tif→mrc
for j in `seq -f %01g 0 12`
do
    tif2mrc mark${j}/t*.tiff mark${j}_mit.mrc
done