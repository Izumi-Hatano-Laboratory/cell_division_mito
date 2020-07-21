#!/bin/sh
i=1
for i in `seq -f %03g 1 400`
do
    echo $i
    convert "p0-${i}.png" "t0-${i}.tiff"
    convert "p1-${i}.png" "t1-${i}.tiff"
    convert "p2-${i}.png" "t2-${i}.tiff"
    convert "p3-${i}.png" "t3-${i}.tiff"
done

tif2mrc t0-*.tiff cluster.mrc
tif2mrc t1-*.tiff small_mit.mrc
tif2mrc t2-*.tiff middle_mit.mrc
tif2mrc t3-*.tiff large_mit.mrc
