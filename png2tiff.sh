#!/bin/sh
i=1
for i in `seq -f %03g 1 400`
do
    echo $i
    convert "mark0/p${i}.png" "mark0/t${i}.tiff"
    convert "mark1/p${i}.png" "mark1/t${i}.tiff"
    convert "mark2/p${i}.png" "mark2/t${i}.tiff"
    convert "mark3/p${i}.png" "mark3/t${i}.tiff"
done

tif2mrc mark0/t*.tiff mark0_mit.mrc
tif2mrc mark1/t*.tiff mark1_mit.mrc
tif2mrc mark2/t*.tiff mark2_mit.mrc
tif2mrc mark3/t*.tiff mark3_mit.mrc
