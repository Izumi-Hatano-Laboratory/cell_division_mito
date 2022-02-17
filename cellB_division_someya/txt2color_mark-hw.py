# -*- coding: utf-8 -*- 

import numpy as np
from PIL import Image, ImageDraw
import os

srcDir = 'D:/cell_division/txt2color'
dirlist = []
filelist1 = []
for root, dirs, files in os.walk(os.path.join(srcDir, 'mark_txt')):
    for file in files:
        filelist1.append(file)
n_txt = len(filelist1)

n256 = 10
for i in range(n256):
    if os.path.isdir("mark{0}".format(i))==False:
        os.mkdir("mark{0}".format(i))

g_max = 0
for j in range(1, n_txt+1):
    mark = np.loadtxt('mark_txt/{:03d}.txt'.format(j), dtype = 'int', delimiter=',')
    mark = mark.reshape(1024,1024)
    print(mark.shape)
    imgsize = (1024,1024)
    img = []
    for i in range(0, n256):
       img.append(Image.new('RGB',imgsize))
    for x in range(imgsize[0]):
        for y in range(imgsize[1]):
            g = mark[y,x]
	    g_max = max([g, g_max])
	    l = (g-1)//255
	    k = (g-1)%255+1
	    for i in range(0,n256):
	     	if i == l:
	    	   img[i].putpixel((x,y),(k,k,k))
		else:
	           img[i].putpixel((x,y),(0,0,0))
    for i in range(0,n256):
           img[i].save('mark{:01d}/p{:03d}.png'.format(i,j), fmt='%s')
print(g_max)
