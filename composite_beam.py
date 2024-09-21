import numpy as np
from numba import jit
@jit()

import numba as nb
import numpy as np
from numba import jit
@jit()

def compbeam(array, hil, samp=100, dx=0.004, sl_low=-2, sl_high=2, sl_ds=0.05, win=4, overlap=0.8):
    
    nch, npts = array.shape
    afinal = np.zeros(npts)
    winlen = win*samp
    sl_2 = np.zeros(2)
    sl_2[0] = np.abs(sl_low)
    sl_2[1] = np.abs(sl_high)
    sl_max = np.max(sl_2)
    n1 = round((0.5*nch*dx*sl_max*samp))
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
                atemp = array[ch,win+dt:win+winlen+dt]
    #             print(win)
    #             print(atemp.size)
                norm = np.sqrt(np.sum(np.square(atemp)))
                atemp = atemp/norm
                cohsum = cohsum + hil[ch,win+dt:win+winlen+dt]
                ssum = ssum + atemp
    #             print(ssum.size)
            ssum = ssum/nch
            cohsum = cohsum/nch 
            ssum = ssum * np.abs(cohsum)**2
#             print(np.sum(ssum**2))
            vesp[k,rc7+5]=np.sum(ssum**2)
        
        
        
        best_idx = np.argmax(np.abs(vesp[:,rc7]))
        
        ssum[:]=0
        cohsum[:]=0;
        ssum2=0
        emat = np.zeros(nch)
        
        for ch, dx in enumerate(x_array):
            dt = int(slowness[best_idx]*dx*samp)
            atemp = array[ch,win+dt:win+winlen+dt]
            ssum2 = ssum2 + np.mean(atemp**2)
            emat[ch] = np.mean(atemp**2)
            norm = np.sqrt(np.sum(np.square(atemp)))
            atemp = atemp/norm
            cohsum = cohsum + hil[ch,win+dt:win+winlen+dt]
            ssum = ssum + atemp
            
        ssum=ssum/nch
        cohsum=cohsum/nch
        ssum = ssum * np.abs(cohsum)**2
        ssum2=np.median(emat)
        ssum=ssum * ssum2
        
        afinal[win:win+winlen]=afinal[win:win+winlen] + ssum
        
    afinal=afinal/5
    

    return vesp, afinal
