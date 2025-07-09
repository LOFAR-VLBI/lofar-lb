import os,glob,numpy as np,matplotlib,h5py
from matplotlib import pyplot as plt
from astropy.time import Time
from astropy.visualization.wcsaxes import SphericalCircle
from astropy.coordinates import SkyCoord
import astropy.io.fits as pyfits; from astropy import units as u
from astropy.wcs import WCS
ANTPLOT = 'DE605'
ANTREF = 'CS001'
#LBCSFILE = 'lbcs_stats.sum'
#LBCSHTTP = 'http://www.jb.man.ac.uk/~njj/lbcs_stats.sum'
POL = 0
CHAN = 50

def movie (title='zoom',fps=4):
#   ---- commented out (gives higher quality but no mencoder available ----
#    a=np.sort(glob.glob(title+'*.png'))
#    for i in range(len(a)):
#        os.system('convert %s %s%06d.jpeg'%(a[i],title,i))
#    command = "mencoder \"mf://%s*.jpeg\" -mf fps=%d -o %s.avi -ovc lavc -lavcopts vcodec=mpeg4:vbitrate=4800" %(title,fps,title)
#    os.system(command)
#   ---- use convert instead ------------
#    os.system('convert '+title+'*.png h5movie.mpg')
#ffmpeg -framerate 4 -pattern_type glob -i 'zoom*.jpg' -c:v libx264 -pix_fmt yuv420p montage.mp4
    pass

def getpcol(phase):
    phase = 3.0*(phase%(2*np.pi))/(2.*np.pi)
    if phase<1.0:
        return (1.0-phase,phase,0.0)
    elif 1.0<=phase<2.0:
        return (0.0,2.0-phase,phase-1.0)
    else:
        return (phase-2.0,0.0,3.0-phase)

def mkfits (rasiz,decsiz,imsiz,pixsiz):
    hdu=pyfits.PrimaryHDU(np.zeros((int(imsiz),int(imsiz))))
    hdu.header.update({'CTYPE1':'RA---SIN'})
    hdu.header.update({'CRVAL1':rasiz})
    hdu.header.update({'CRPIX1':imsiz/2.})
    hdu.header.update({'CDELT1':-pixsiz})
    hdu.header.update({'CTYPE2':'DEC--SIN'})
    hdu.header.update({'CRVAL2':decsiz})
    hdu.header.update({'CRPIX2':imsiz/2.})
    hdu.header.update({'CDELT2':pixsiz})
    hdu.header.update({'EQUINOX':2000.0})
    hdu.header.update({'EPOCH':2000.0})
    os.system('rm temp.fits')
    hdu.writeto('temp.fits')

def get_h5list ():   # adjust for the system used
    os.system('find * -print |grep merged|grep h5 >temp')
    h5files = np.loadtxt('temp',dtype='str')
    return np.sort(h5files)

h5list = get_h5list()
#try:
#    lbcs = np.loadtxt(LBCSFILE,dtype='str')
#except:
#    os.system('wget '+LBCSHTTP)
#    lbcs = np.loadtxt(LBCSFILE,dtype='str')

pos=np.array([]);phase=np.array([]);ut=np.array([])
for h5 in h5list:
    h = h5py.File(h5)
    newpos = np.asarray(h['sol000/source'][0][1],dtype='float')
    ant = np.array(h['sol000/antenna'])
    for i in range(len(ant)):
        if str(ant[i][0],'utf8')[:5] == ANTPLOT:
            idxant = i
        if str(ant[i][0],'utf8')[:5] == ANTREF:
            idxref = i
    newphase = np.asarray(h['sol000/phase000/val'],dtype='float')[:,CHAN,idxant,POL,0]
    newref = np.asarray(h['sol000/phase000/val'],dtype='float')[:,CHAN,idxref,POL,0]
    newut = np.asarray(h['sol000/phase000']['time'],dtype='float')/86400.0
    try:
        ut = np.vstack((ut,newut))
        pos = np.vstack((pos,newpos))
        phase = np.vstack((phase,newphase-newref))
    except:
        ut = np.copy(newut)
        pos = np.copy(newpos)
        phase = np.copy(newphase-newref)

REFSRC = len(phase)//2
for i in range(len(phase)):
    if i!=REFSRC:
        phase[i]-=phase[REFSRC]

phase[REFSRC]*=0.0
np.putmask(pos[:,0],pos[:,0]<0.0,pos[:,0]+2*np.pi)
pos *= 180./np.pi
rawidth = pos[:,0].max()-pos[:,0].min()
decwidth = pos[:,1].max()-pos[:,1].min()
ra = (pos[:,0].max()+pos[:,0].min())/2.
dec = (pos[:,1].max()+pos[:,1].min())/2.
width = 1.5*max(decwidth,rawidth*np.cos(dec))
mkfits (ra,dec,600,width/600.0)
#for nframe in range(phase.shape[1]):
os.system('rm zoom*.png;rm zoom*jpeg')
print(phase.shape)
for nframe in range(min(phase.shape[1],500)):
    print('Frame %d'%nframe)
    hdu = pyfits.open('temp.fits')[0]
    ax = plt.subplot(1,1,1,xticks=[],yticks=[],projection=WCS(hdu.header))
    plt.imshow(pyfits.getdata('temp.fits'),cmap=matplotlib.cm.gray_r)
#    ax.coords[0].set_axislabel(False)
#    ax.coords[1].set_axislabel(False)
#    ax.coords[0].set_ticks_visible(False)
#    ax.coords[1].set_ticklabel_visible(False)
#    ax.coords[0].set_ticklabel_visible(False)
#    ax.coords[1].set_ticks_visible(False)
    for i in range(len(phase)):
        ax.scatter(pos[i,0],pos[i,1],s=100,facecolor=getpcol(phase[i,nframe]),transform=ax.get_transform('fk5'))
#        thiscoord = SkyCoord(pos[i,0]*u.deg,pos[i,1]*u.deg)
#        print (thiscoord,phase[i,nframe],getpcol(phase[i,nframe]))
#        ell = SphericalCircle (thiscoord,0.01*u.deg,\
#                               facecolor = getpcol(phase[i,nframe]))
#    ax.add_artist(ell)
    plt.title(ANTPLOT)
    plt.savefig('zoom%06d.png'%nframe)
    plt.clf()

#    gc = aplpy.FITSFigure('temp.fits')
#    gc.add_grid()
#    gc.grid.set_color('black')
#    for i in range(len(phase)):
#        gc.show_circles(pos[i,0],pos[i,1],0.1,color='k',\
#                    facecolor=getpcol(phase[i,nframe]))
#    gc.set_title(Time(ut[0,nframe],format='mjd').fits)
#    gc.save('zoom%06d.png'%nframe)
#    gc.close()

#movie (title='zoom',fps=10)
