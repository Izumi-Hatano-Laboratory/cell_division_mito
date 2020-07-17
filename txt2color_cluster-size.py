import numpy as np
from PIL import Image, ImageDraw

g = [0,0,0,0]
for j in range(1,401):
    cluster = np.loadtxt('cluster/{:03d}.txt'.format(j), dtype = 'int', delimiter=',')
    cluster = cluster.reshape(1024,1024)
    size1 = np.loadtxt('size_0_100000/{:03d}.txt'.format(j), dtype = 'int', delimiter=',')
    size1 = size1.reshape(1024,1024)
    size2 = np.loadtxt('size_100000_200000/{:03d}.txt'.format(j), dtype = 'int', delimiter=',')
    size2 = size2.reshape(1024,1024)
    size3 = np.loadtxt('size_200000_upper/{:03d}.txt'.format(j), dtype = 'int', delimiter=',')
    size3 = size3.reshape(1024,1024)
#    print(cluster.shape)
    imgsize = (1024,1024)

    img = []
    nimg = 4
    for i in range(0, nimg):
       img.append(Image.new('RGB',imgsize))
    for x in range(imgsize[0]):
        for y in range(imgsize[1]):
            g[0] = cluster[y,x]%255
            g[1] = size1[x,y]%255
            g[2] = size2[x,y]%255
            g[3] = size3[x,y]%255
	    for i in range(0,nimg):
	    	img[i].putpixel((x,y),(g[i],g[i],g[i]))
    for i in range(0,nimg):
           img[i].save('cluster/p{:01d}-{:03d}.png'.format(i,j), fmt='%s')


