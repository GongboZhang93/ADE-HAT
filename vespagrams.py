import numpy as np
from numba import jit
@jit()

def vespagram(array, hil, samp=100, dx=0.004, sl_low=-2, sl_high=2, sl_ds=0.05, win=4, overlap=0.8):
    nch, npts = data.shape
    winlen = win*samp
    n1 = round((0.5*nch*dx*sl_high*samp))
    n2 = n1+winlen
    interval = round((1-0.8)*400)
    slowness = np.arange(sl_low, sl_high+sl_ds, sl_ds)
    nslow = slowness.size
    win_idx = np.arange(n1,npts-n2,interval)
    nwin = win_idx.size
    nwin0 = np.arange(0,npts,interval).size
    x_array = np.arange(nch) * dx
    x_array -= x_array.mean()  #km
#     print(data.shape)
#     print(data.shape)
    
    vesp = np.zeros((nslow, nwin0))
    for rc7, win in enumerate(win_idx):
        print(rc7)
        for k, s, in enumerate(slowness):
            ssum = np.zeros(winlen)
            cohsum = np.zeros(winlen, dtype = nb.complex128)
            for ch, dx in enumerate(x_array):
                dt = int(s*dx*samp)
                atemp = data[ch,win+dt:win+winlen+dt]
    #             print(win)
    #             print(atemp.size)
                norm = np.sqrt(np.sum(np.square(atemp)))
                atemp = atemp/norm
                cohsum = cohsum + data_hil[ch,win+dt:win+winlen+dt]
                ssum = ssum + atemp
    #             print(ssum.size)
            ssum = ssum/nch
            cohsum = cohsum/nch 
            ssum = ssum * np.abs(cohsum)**2
#             print(np.sum(ssum**2))
            vesp[k,rc7+5]=np.sum(ssum**2)
        
    return vesp

