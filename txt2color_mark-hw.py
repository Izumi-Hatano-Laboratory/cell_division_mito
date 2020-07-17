import numpy as np
from PIL import Image, ImageDraw

g_max = 0
for j in range(1,401):
    mark = np.loadtxt('mark_txt/mark-hw{:03d}.txt'.format(j), dtype = 'int', delimiter=',')
    mark = mark.reshape(1024,1024)
    print(mark.shape)
    imgsize = (1024,1024)

    img = []
    n256 = 2
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
#       if j < 10:
#           img[i].save('mark'+str(i)+'/p00'+str(j)+'.png', fmt='%s')
#       elif j < 100:
#           img[i].save('mark'+str(i)+'/p0'+str(j)+'.png', fmt='%s')
#       elif j < 1000:
#           img[i].save('mark'+str(i)+'/p'+str(j)+'.png', fmt='%s')
print(g_max)

