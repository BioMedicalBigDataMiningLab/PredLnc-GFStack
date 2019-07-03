
import numpy as np
#SNR
def SNR(seq):
    Alist = DFT_nlist(seq,'A')
    Glist = DFT_nlist(seq,'G')
    Clist = DFT_nlist(seq,'C')
    Tlist = DFT_nlist(seq,'T')
    Plist = [abs(Alist[i])**2+abs(Glist[i])**2+abs(Clist[i])**2+abs(Tlist[i])**2 for i in range(len(seq))]
    R = Plist[round(len(seq)/3)]/(sum(Plist)/len(seq))
    return R

def DFT_nlist(seq,basic):
    nlist = []
    for i in range(len(seq)):
        if seq[i]==basic:
            nlist.append(1)
        else:
            nlist.append(0)
    nlist = np.fft.fft(nlist)
    return nlist
