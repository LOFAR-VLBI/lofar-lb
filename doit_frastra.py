import os,sys,numpy as np
#os.system('source /home/njj/.cshrc')
#os.system('source /aips/LOGIN.CSH')
#os.system('sed "s/source[i]/sources[i]/g" proc_pmkmap.py >proc_pmkmap1.py')
#os.system('grep sources proc_pmkmap1.py')
#os.system('which parseltongue')
#os.system('parseltongue proc_pmkmap1.py ../newdata/44059320F.fits')
os.system('rm labdem.pdf')
os.system('wget https://www.dropbox.com/scl/fi/8ncyaeecr4j9noz8i53xh/labdem.pdf')
os.system('mv -f labdem.pdf /home/njj/public_html/labdem.pdf')
os.system('rm doit_frastra.py')
