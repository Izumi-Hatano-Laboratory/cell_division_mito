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
while [ $i -lt  766 ] #255*3+1  �Ō�̓~�g�R���h���A��
do
    echo $i
    j=`expr $i - 510`
    imodauto "-E" "$j" "mark2_mit.mrc" "mark_tmp.mod"
    imodjoin "-c" "mark_mit.mod" "mark_tmp.mod" "mark_mit.mod"
    i=`expr $i + 1`
done
while [ $i -lt  1021 ] #255*4+1  �Ō�̓~�g�R���h���A��
do
    echo $i
    j=`expr $i - 765`
    imodauto "-E" "$j" "mark3_mit.mrc" "mark_tmp.mod"
    imodjoin "-c" "mark_mit.mod" "mark_tmp.mod" "mark_mit.mod"
    i=`expr $i + 1`
done
while [ $i -lt  1276 ] #255*5+1  �Ō�̓~�g�R���h���A��
do
    echo $i
    j=`expr $i - 1020`
    imodauto "-E" "$j" "mark4_mit.mrc" "mark_tmp.mod"
    imodjoin "-c" "mark_mit.mod" "mark_tmp.mod" "mark_mit.mod"
    i=`expr $i + 1`
done
while [ $i -lt  1531 ] #255*6+1  �Ō�̓~�g�R���h���A��
do
    echo $i
    j=`expr $i - 1275`
    imodauto "-E" "$j" "mark5_mit.mrc" "mark_tmp.mod"
    imodjoin "-c" "mark_mit.mod" "mark_tmp.mod" "mark_mit.mod"
    i=`expr $i + 1`
done
while [ $i -lt  1786 ] #255*7+1  �Ō�̓~�g�R���h���A��
do
    echo $i
    j=`expr $i - 1530`
    imodauto "-E" "$j" "mark6_mit.mrc" "mark_tmp.mod"
    imodjoin "-c" "mark_mit.mod" "mark_tmp.mod" "mark_mit.mod"
    i=`expr $i + 1`
done
while [ $i -lt  2041 ] #255*8+1  �Ō�̓~�g�R���h���A��
do
    echo $i
    j=`expr $i - 1785`
    imodauto "-E" "$j" "mark7_mit.mrc" "mark_tmp.mod"
    imodjoin "-c" "mark_mit.mod" "mark_tmp.mod" "mark_mit.mod"
    i=`expr $i + 1`
done
while [ $i -lt  2296 ] #255*9+1  �Ō�̓~�g�R���h���A��
do
    echo $i
    j=`expr $i - 2040`
    imodauto "-E" "$j" "mark8_mit.mrc" "mark_tmp.mod"
    imodjoin "-c" "mark_mit.mod" "mark_tmp.mod" "mark_mit.mod"
    i=`expr $i + 1`
done
while [ $i -lt  2551 ] #255*10+1  �Ō�̓~�g�R���h���A��
do
    echo $i
    j=`expr $i - 2295`
    imodauto "-E" "$j" "mark9_mit.mrc" "mark_tmp.mod"
    imodjoin "-c" "mark_mit.mod" "mark_tmp.mod" "mark_mit.mod"
    i=`expr $i + 1`
done
while [ $i -lt  2806 ] #255*11+1  �Ō�̓~�g�R���h���A��
do
    echo $i
    j=`expr $i - 2550`
    imodauto "-E" "$j" "mark10_mit.mrc" "mark_tmp.mod"
    imodjoin "-c" "mark_mit.mod" "mark_tmp.mod" "mark_mit.mod"
    i=`expr $i + 1`
done
while [ $i -lt  3061 ] #255*12+1  �Ō�̓~�g�R���h���A��
do
    echo $i
    j=`expr $i - 2805`
    imodauto "-E" "$j" "mark11_mit.mrc" "mark_tmp.mod"
    imodjoin "-c" "mark_mit.mod" "mark_tmp.mod" "mark_mit.mod"
    i=`expr $i + 1`
done
while [ $i -lt  3316 ] #255*13+1  �Ō�̓~�g�R���h���A��
do
    echo $i
    j=`expr $i - 3060`
    imodauto "-E" "$j" "mark12_mit.mrc" "mark_tmp.mod"
    imodjoin "-c" "mark_mit.mod" "mark_tmp.mod" "mark_mit.mod"
    i=`expr $i + 1`
done
<< COMMENTOUT
while [ $i -lt  3571 ] #255*14+1  �Ō�̓~�g�R���h���A��
do
    echo $i
    j=`expr $i - 3315`
    imodauto "-E" "$j" "mark13_mit.mrc" "mark_tmp.mod"
    imodjoin "-c" "mark_mit.mod" "mark_tmp.mod" "mark_mit.mod"
    i=`expr $i + 1`
done
while [ $i -lt  3826 ] #255*15+1  �Ō�̓~�g�R���h���A��
do
    echo $i
    j=`expr $i - 3570`
    imodauto "-E" "$j" "mark14_mit.mrc" "mark_tmp.mod"
    imodjoin "-c" "mark_mit.mod" "mark_tmp.mod" "mark_mit.mod"
    i=`expr $i + 1`
done
while [ $i -lt  4081 ] #255*16+1  �Ō�̓~�g�R���h���A��
do
    echo $i
    j=`expr $i - 3825`
    imodauto "-E" "$j" "mark15_mit.mrc" "mark_tmp.mod"
    imodjoin "-c" "mark_mit.mod" "mark_tmp.mod" "mark_mit.mod"
    i=`expr $i + 1`
done
COMMENTOUT
#imodmesh -c mark_mit.mod
imodinfo -c -f mark_mit_info.txt mark_mit.mod

#imodauto:mrcからモデルを作成
#imodjoin:複数のimodモデルをくっつける
#imodmesh:modからメッシュ切る
#imodinfo:情報を貼る