import numpy as np
from PIL import Image, ImageDraw

for i in range(1,401):
    img   = Image.open('./ensembled_full2000_InsideLearned/Segmented_{:04d}.png'.format(i))
    img_s = Image.open('./ensembled_full20000_ContourLearned/Segmented_{:04d}.png'.format(i))
    img = img.convert('RGB')
    img_s = img_s.convert('RGB')
    
    width, height = img_s.size
    
    contour = np.array([[img_s.getpixel((y,x)) for y in range(height)] for x in range(width)])
    contour = np.ravel(contour[0:1024,0:1024,1])
    mit = np.array([[img.getpixel((y,x)) for y in range(height)] for x in range(width)])
    mit = np.ravel(mit[0:1024,0:1024,1])

    np.savetxt('out_txt/{:03d}.txt'.format(i), contour, fmt='%s', delimiter=',')
    np.savetxt('in_txt/{:03d}.txt'.format(i), mit, fmt='%s', delimiter=',')
