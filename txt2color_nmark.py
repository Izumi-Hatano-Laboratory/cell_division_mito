import numpy as np
from PIL import Image, ImageDraw
#lists = [1, 21, 41]
#for j in lists:
for j in range(1,11):
    nmark = np.loadtxt('nmark/{:03d}.txt'.format(j), dtype = 'int', delimiter=',')
    nmark = nmark.reshape(1024,1024)
    print(nmark.shape)
    print(nmark.dtype)
    print(nmark.ndim)
    
   
    imgsize = (1024,1024)#nmark.size#[1],nmark.size[0])
    img = Image.new('RGB',imgsize)
    for x in range(imgsize[0]):
        for y in range(imgsize[1]):
            g = nmark[y,x]
	    r = g//255
	    b = g//(255*255)
	    g = g%255
            img.putpixel((x,y),(g,g,g))

    img.save('nmark/{:03d}.png'.format(j), fmt='%s')

