Folder Structure
./in_txt/???.txt   #ミトコンドリア内部を学習した結果を連番txtにしたものを格納
./out_txt/???.txt  #ミトコンドリア輪郭を学習した結果を連番txtにしたものを格納
./nucleus/         #中間結果用フォルダを用意しておく　抽出した核
./nmark/           #　　　　　　　　　　　　　　　　　核毎に別色付け
./smark/           #　　　　　　　　　　　　　　　　　輪郭近くまでの色付け
./mark/            #　　　　　　　　　　　　　　　　　ミトコンドリア色づけ
./mark_txt/        #結果フォルダ　小さすぎるミトコンドリア等を除いた結果格納
./mark0/           #png保存用フォルダ．下記n256数に合わせて用意すること．
./mark1/
./mark2/


実行プロセス
python color2txt.py  # in_txt, out_txt を準備する．
make                 # marker_3dout.f をコンパイルして実行ファイルmarkerを作成
./marker             # mark_txt に結果が入る．1hほどかかる
python txt2color.py  # txtをpngに変える．
                     # ミトコンドリアの数が255を超えると同じ色が出てしまうので
                     # それを分類するために255おきに別フォルダに入れる
                     # そのためのフォルダがmark0, mark1, mark2 ... 
                     # ミトコンドリア数/255 +1 = n256 とすること．
png2tiff.sh          # pngをtiffにして，mrcファイル（imod用)を作成
                     # 3dmod インストールされたマシンで実行すること．
                     # 上記フォルダ数（n256)に合わせてスクリプト変更のこと．
sh imodauto_auto.sh  # imodauto -E 値 file.mrc tmp.mod
                     #   mrcより値の領域をモデルとして保存．境界を自動作成
                     # imodjoin -c tmp.mod mit.mod mit.mod
                     #   tmp.mod と mit.mod をあわせたものをmit.mod として上書
                     #   -cは毎回違う色でという指定．
