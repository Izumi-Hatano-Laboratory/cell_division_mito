Folder Structure
./in_txt/???.txt   #�~�g�R���h���A�������w�K�������ʂ�A��txt�ɂ������̂��i�[
./out_txt/???.txt  #�~�g�R���h���A�֊s���w�K�������ʂ�A��txt�ɂ������̂��i�[
./nucleus/         #���Ԍ��ʗp�t�H���_��p�ӂ��Ă����@���o�����j
./nmark/           #�@�@�@�@�@�@�@�@�@�@�@�@�@�@�@�@�@�j���ɕʐF�t��
./smark/           #�@�@�@�@�@�@�@�@�@�@�@�@�@�@�@�@�@�֊s�߂��܂ł̐F�t��
./mark/            #�@�@�@�@�@�@�@�@�@�@�@�@�@�@�@�@�@�~�g�R���h���A�F�Â�
./mark_txt/        #���ʃt�H���_�@����������~�g�R���h���A�������������ʊi�[
./mark0/           #png�ۑ��p�t�H���_�D���Ln256���ɍ��킹�ėp�ӂ��邱�ƁD
./mark1/
./mark2/


���s�v���Z�X
python color2txt.py  # in_txt, out_txt ����������D
make                 # marker_3dout.f ���R���p�C�����Ď��s�t�@�C��marker���쐬
./marker             # mark_txt �Ɍ��ʂ�����D1h�قǂ�����
python txt2color.py  # txt��png�ɕς���D
                     # �~�g�R���h���A�̐���255�𒴂���Ɠ����F���o�Ă��܂��̂�
                     # ����𕪗ނ��邽�߂�255�����ɕʃt�H���_�ɓ����
                     # ���̂��߂̃t�H���_��mark0, mark1, mark2 ... 
                     # �~�g�R���h���A��/255 +1 = n256 �Ƃ��邱�ƁD
png2tiff.sh          # png��tiff�ɂ��āCmrc�t�@�C���iimod�p)���쐬
                     # 3dmod �C���X�g�[�����ꂽ�}�V���Ŏ��s���邱�ƁD
                     # ��L�t�H���_���in256)�ɍ��킹�ăX�N���v�g�ύX�̂��ƁD
sh imodauto_auto.sh  # imodauto -E �l file.mrc tmp.mod
                     #   mrc���l�̗̈�����f���Ƃ��ĕۑ��D���E�������쐬
                     # imodjoin -c tmp.mod mit.mod mit.mod
                     #   tmp.mod �� mit.mod �����킹�����̂�mit.mod �Ƃ��ď㏑
                     #   -c�͖���Ⴄ�F�łƂ����w��D
