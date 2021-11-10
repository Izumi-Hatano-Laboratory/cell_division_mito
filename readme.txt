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

～～染谷加筆～～
＜color2txt.py＞
深層学習結果をテキストデータin_txt, out_txtに変換
srcDirを適宜変更
＊フォルダ構造＊
srcDir----in		    小器官内部の画像、深層学習結果のグレースケールそのまま
	      |-in_txt		画像のtxtデータ
	      |-out		    小器官の輪郭画像
	      |-out_txt		画像のtxtデータ

＜marker_contourOut13.f＞
分割アルゴリズム。mark_txtに結果が入る。
画像の枚数depth(l.15)、適宜変更
二値化のしきい値(l.172~182)、ミトコンドリア内の内膜が消えるしきい値と、差分をとったあとに個々を分離するしきい値にそれぞれ適当な値を設定。

＜txt2color_mark-hw.py＞
分割後のtxtデータをpngにする。時間かかる。256個ずつしかミトコンドリアを同じ画像に表示できないため、出力フォルダを(ミトコンドリアの数÷256)個作る
srcDir(l.7), n256(>=ミトコンドリア数//255 +1)(l.15)、を適宜変更
＊フォルダ構造＊
srcDir-----mark_txt	markerで分割した結果
	       |-mark0		.pngの出力場所、n256個分
	       |-mark1		
         |-　・		
         |-　・		
         |-mark9		(n256=10のとき)
	      
＜png2tiff.sh＞
pngを更にtiff, mrcに変換
3dmodのインストールされたマシンで行う
画像枚数(l.3)、markフォルダの数(l.8,15)、変更

＜imod_auto_auto.sh＞
mark**_mit.modをすべて結合する
.modファイルの数だけwhile do ~ doneを繰り返す。
