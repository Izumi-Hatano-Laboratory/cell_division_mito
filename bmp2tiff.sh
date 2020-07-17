#!/bin/sh
i=1
for i in `seq -f %03g 1 400`
do
    echo $i
    convert "${i}.bmp" "${i}.tiff"
done
