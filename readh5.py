import numpy as np,os,sys,h5py,scipy,glob
import matplotlib; from matplotlib import pyplot as plt
from scipy.optimize import curve_fit
import tables
matplotlib.rcParams['image.origin'] = 'lower'
matplotlib.rcParams['image.interpolation'] = 'nearest'

def sort_by_number(antennas):
    iant = []
    for i in antennas:
        nidx = int(i[2:5])
        if i[:2]=='CS':
            nidx-=2000
        if i[:2]=='RS':
            nidx-=1000
        iant.append(nidx)
    return antennas[np.argsort(iant)]

def metric_test (a):
    np.putmask(a,np.isnan(a),-np.pi*2*np.random.random())
    return np.median(np.abs(np.gradient(np.unwrap(a))))

def plotant(data,antn,title=None,t1=0,t2=-1):
    if t2==-1:
        t2=data.shape[0]
    aspect=50/(t2-t1)
#    if len(data.shape)==4:
#        data = data[:,:,:,0]-data[:,:,:,1]
    data = data[:,:,:,0]
    plt.clf()
    plt.imshow(data[t1:t2,antn,:],cmap=matplotlib.cm.jet,aspect=aspect,\
               vmin=-np.pi,vmax=np.pi)
    if title is not None:
        plt.title(title)
#    plt.colorbar()
    plt.savefig(str(title)+'.png',bbox_inches='tight')
#    plt.show()

def h5_read (filename):
    h5 = tables.open_file(filename)
    try:
        h5_d=h5.root['target']['TGSSscalarphase_final']['val'].read()
        h5_a=h5.root['target']['TGSSscalarphase_final']['ant'].read()
        h5_a=np.asarray(h5_a,dtype='str')
        h5_f=h5.root['target']['TGSSscalarphase_final']['freq'].read()
        h5_t=h5.root['target']['TGSSscalarphase_final']['time'].read()
        h5_type = 'scalar phase'
    except:
        h5_d=h5.root['target']['TGSSphase']['val'].read()
        h5_a=h5.root['target']['TGSSphase']['ant'].read()
        h5_a=np.asarray(h5_a,dtype='str')
        h5_f=h5.root['target']['TGSSphase']['freq'].read()
        h5_t=h5.root['target']['TGSSphase']['time'].read()
        h5_type = 'phase'
    sys.stdout.write ('*** %s ***\nType %s: %d channels, %d antennas, int. time %.3fs, total %.3fhr\n'\
           %(filename,h5_type,h5_d.shape[2],h5_d.shape[1],\
             h5_t[1]-h5_t[0],(h5_t[1]-h5_t[0])*h5_d.shape[0]/3600))
    return h5_d,h5_a,h5_f,h5_t[1]-h5_t[0]

def proc_solutions (filename,chan='middle',verbose=True,anttoplot='RS310HBA',\
                    doplot=False):
    if not 'h5' in filename:
        print ('Not an h5 file, exiting')
        return
    data, antennas, freqs, inttime = h5_read(filename)
    if doplot:
        print(antennas[0],np.argwhere(antennas==anttoplot)[0][0])
        plotant (data, np.argwhere(antennas==anttoplot)[0][0],\
                 title=filename.split('.')[0]+anttoplot)
    if chan=='middle':
        chan=int(len(freqs)/2)
    antennas_n = sort_by_number(antennas)
    pic = np.zeros((len(antennas),8),dtype='float')
    thr = [0.13,0.2]
    for i, antenna in enumerate(antennas_n):
        antenna_idx = np.argwhere(antennas==antenna)[0][0]
        if verbose:
            sys.stdout.write ('%8s '%antenna[:5])
        for j,tstart in enumerate(np.arange(0.0, 480.0, 60.0)):
            tend = tstart + 60.
            istart,iend = int(tstart*60/inttime), int(tend*60/inttime)
            iend = min(data.shape[0], iend)
            if istart>data.shape[0] or iend-istart<10:
                if verbose:
                    sys.stdout.write ('-----')
                continue
            if len(data.shape)==4:
                val = metric_test (data[istart:iend,antenna_idx,chan,0])
            else:
                val = metric_test (data[istart:iend,antenna_idx,chan])
            pic[i,j] = val
            if verbose:
                sys.stdout.write ('%.3f '%(val))
        if verbose:
            sys.stdout.write('\n')

    plt.clf()
    plt.imshow(pic,cmap=matplotlib.cm.jet,aspect='auto',vmin=0,vmax=0.35)
    plt.colorbar()
    plt.title(filename.split('.')[0]+'chan%d'%chan)
    plt.savefig(filename.split('.')[0]+'chan%d'%chan,bbox_inches='tight')
    return antennas_n, pic

def action (ant, a, threshold=0.15, frac=0.33):
    remove_ants = []
    for i in range(len(ant)):
        gvals = a[i][a[i]>0.0]
        if len(gvals) and len(gvals[gvals>threshold])/len(gvals)>frac:
            remove_ants.append(ant[i])
    print('Recommended deletions: ',remove_ants)
    

if len(sys.argv)<2:
    for i in glob.glob('L*solutions.h5'):
        ant, pic = proc_solutions(i,verbose=False,doplot=True)
        action (ant, pic)
else:
     ant, pic = proc_solutions('L%s_solutions.h5'%sys.argv[1],\
                               verbose=False,doplot=True)
     action (ant, pic)


# function limbo

def phase_fn (freq, a, b, c):
    return a+b*freq+c/freq   # constant, clock term, disp term

def fit_freq (freq, data):
    popt,pcov = curve_fit(phase_fn,freq,data)
    return popt,phase_fn (freq, popt[0],popt[1],popt[2])

def print_dataset_info(name, obj, nval=1):
    """Recursively prints details of datasets, including nested ones."""
    if isinstance(obj, h5py.Dataset):  # If it's a dataset
        print(f"Dataset '{name}': Shape = {obj.shape}, Type = {obj.dtype}")
        print(f"First few values: {obj[:nval]}")
    elif isinstance(obj, h5py.Group):  # If it's a group, explore further
        print(f"Group '{name}': (Contains nested datasets or subgroups)")
        for key in obj.keys():
            print_dataset_info(f"{name}/{key}", obj[key])  # Recursive call

def print_info(filename,field='target',nval=1):
    with h5py.File(filename,'r') as f:
        print_dataset_info(field,f[field],nval=nval)


def print_dataset_info(name, obj, nval=1):
    """Recursively prints details of datasets, including nested ones."""
    if isinstance(obj, h5py.Dataset):  # If it's a dataset
        print(f"Dataset '{name}': Shape = {obj.shape}, Type = {obj.dtype}")
        print(f"First few values: {obj[:nval]}")
    elif isinstance(obj, h5py.Group):  # If it's a group, explore further
        print(f"Group '{name}': (Contains nested datasets or subgroups)")
        for key in obj.keys():
            print_dataset_info(f"{name}/{key}", obj[key])  # Recursive call

def print_info(filename,field='target',nval=1):
    with h5py.File(filename,'r') as f:
        print_dataset_info(field,f[field],nval=nval)



#data = np.array([[1.2, 2.2,   'o',   's', 2.5],
#                 [1.7,   's', 2.4, 2.9, 1.7],
#                 ['o', 0.9, 0.1, 'NaN', 0.4]])
#data[data == 'NaN'] = -4000           
#data[data == 'o'] = -5000           
#data[data == 's'] = -6000
#data = data.astype(np.float64)

#norm = mpl.colors.Normalize(vmin=0, vmax=3)
#cmap = plt.get_cmap('viridis')

# convert from 0-1:
#datan = norm(data)
# convert to rgba:
#rgba = cmap(datan)
# Fill in colours for the out of range data:
#rgba[data==-4000, :] = [1, 1, 1, 1]
#rgba[data==-5000, :] = [1, 0, 0, 1]
#rgba[data==-6000, :] = [0, 1, 1, 1]

# plot:
#fig, ax = plt.subplots()
#ax.imshow(rgba)
#plt.show()

#enter image description here
