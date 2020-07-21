import numpy as np
from PIL import Image, ImageDraw
im = np.zeros(1024*1024*3, dtype=np.uint8)
im = im.reshape(1024,1024,3)
for j in range(1,401):
    print(j)
    nmark = np.loadtxt('nmark/{:03d}.txt'.format(j), dtype = 'int', delimiter=',')
    nmark = nmark.reshape(1024,1024)
    nmark.transpose()
    nmark = nmark%255
    im = nmark.astype('uint8')
    img = Image.fromarray(im)
    img.save('nmark/{:03d}.png'.format(j), fmt='%s')
