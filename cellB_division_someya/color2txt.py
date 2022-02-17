# coding: UTF-8

import numpy as np
from PIL import Image, ImageDraw
import os

srcDir = 'D:/cell_division/for_color2txt'
dirlist = []
for f in os.listdir(srcDir):
    if os.path.isdir(os.path.join(srcDir, f)):
        dirlist.append(f)

filelist_in = []   
filelist_out = []   
for root, dirs, files in os.walk(os.path.join(srcDir, dirlist[0])):
    for file in files:
        filelist_in.append(file)
for root, dirs, files in os.walk(os.path.join(srcDir, dirlist[2])):
    for file in files:
        filelist_out.append(file)

for i in range(0,len(filelist_in)):
    path = os.path.join(srcDir, dirlist[0], filelist_in[i])
    img = Image.open(path) 
    path = os.path.join(srcDir, dirlist[2], filelist_out[i])
    img_s = Image.open(path)
    img = img.convert('RGB')
    img_s = img_s.convert('RGB')
    width, height = img_s.size
    contour = np.array([[img_s.getpixel((y,x)) for y in range(height)] for x in range(width)])
    contour = np.ravel(contour[0:1024,0:1024,1])  #contourを一次元配列化
    mit = np.array([[img.getpixel((y,x)) for y in range(height)] for x in range(width)])
    mit = np.ravel(mit[0:1024,0:1024,1])
    np.savetxt('out_txt/{:03d}.txt'.format(i), contour, fmt='%s', delimiter=',')
    np.savetxt('in_txt/{:03d}.txt'.format(i), mit, fmt='%s', delimiter=',')
