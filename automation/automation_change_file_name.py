import os

path = os.path.dirname(os.path.realpath(__file__))

files = os.listdir(path)
for xvgfile in files:
    if xvgfile.endswith(".xvg"):
        #print(xvgfile)
        os.rename(os.path.join(path, xvgfile), os.path.join(path, xvgfile[:-4]))

files = os.listdir(path)
for txtfile in files:
    if txtfile[-4:] == '.txt':
        #print(path + xvgfile)
        finn = open(txtfile, 'r')
        a = finn.readlines()
        fout = open(txtfile, 'w')
        b = ''.join(a[13:-1])
        fout.write(b)
