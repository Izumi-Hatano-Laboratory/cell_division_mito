#!/bin/sh
i=1
imodauto "-E" 1 "mark0_mit.mrc" "mark_mit.mod"
i=2
while [ $i -lt 256 ]
do
    echo $i
    imodauto "-E" "$i" "mark0_mit.mrc" "mark_tmp.mod"
    imodjoin "-c" "mark_mit.mod" "mark_tmp.mod" "mark_mit.mod"
    i=`expr $i + 1`
done
while [ $i -lt 511 ] #255*2+1
do
    echo $i
    j=`expr $i - 255`
    imodauto "-E" "$j" "mark1_mit.mrc" "mark_tmp.mod"
    imodjoin "-c" "mark_mit.mod" "mark_tmp.mod" "mark_mit.mod"
    i=`expr $i + 1`
done
while [ $i -lt  766 ] #255*3+1  最後はミトコンドリア数
do
    echo $i
    j=`expr $i - 510`
    imodauto "-E" "$j" "mark2_mit.mrc" "mark_tmp.mod"
    imodjoin "-c" "mark_mit.mod" "mark_tmp.mod" "mark_mit.mod"
    i=`expr $i + 1`
done
while [ $i -lt  1021 ] #255*4+1  最後はミトコンドリア数
do
    echo $i
    j=`expr $i - 765`
    imodauto "-E" "$j" "mark3_mit.mrc" "mark_tmp.mod"
    imodjoin "-c" "mark_mit.mod" "mark_tmp.mod" "mark_mit.mod"
    i=`expr $i + 1`
done

imodmesh -c mark_mit.mod
imodinfo -c -f mark_mit_info.txt mark_mit.mod
