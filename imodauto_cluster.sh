#!/bin/sh
i=1
imodauto "-E" 1 "cluster.mrc" "cluster.mod"
i=2
while [ $i -lt 256 ]
do
    echo $i
    imodauto "-E" "$i" "cluster.mrc" "tmp.mod"
    imodjoin "-c" "cluster.mod" "tmp.mod" "cluster.mod"
    i=`expr $i + 1`
done
imodmesh "-c" "cluster.mod"
imodinfo "-c" "-f" "cluster_info.txt" "cluster.mod"
#SMALL
i=1
imodauto "-E" 1 "small.mrc" "small.mod"
i=2
while [ $i -lt 256 ] #255*2+1
do
    echo $i
    imodauto "-E" "$i" "small.mrc" "tmp.mod"
    imodjoin "-c" "small.mod" "tmp.mod" "small.mod"
    i=`expr $i + 1`
done
imodmesh "-c" "small.mod"
imodinfo "-c" "-f" "small_info.txt" "small.mod"
#MIDDLE
i=1
imodauto "-E" 1 "middle.mrc" "middle.mod"
i=2
while [ $i -lt 256 ] #255*2+1
do
    echo $i
    imodauto "-E" "$i" "middle.mrc" "tmp.mod"
    imodjoin "-c" "middle.mod" "tmp.mod" "middle.mod"
    i=`expr $i + 1`
done
imodmesh "-c" "middle.mod"
imodinfo "-c" "-f" "middle_info.txt" "middle.mod"
#LARGE
i=1
imodauto "-E" 1 "large.mrc" "large.mod"
i=2
while [ $i -lt 256 ] #255*2+1
do
    echo $i
    imodauto "-E" "$i" "large.mrc" "tmp.mod"
    imodjoin "-c" "large.mod" "tmp.mod" "large.mod"
    i=`expr $i + 1`
done
imodmesh "-c" "large.mod"
imodinfo "-c" "-f" "large_info.txt" "large.mod"
