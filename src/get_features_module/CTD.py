#Return dictionary
def CTD(seq):
    seq = seq.upper()
    n = len(seq)
    names = ['A','G','C','T','AA', 'AAA', 'AAC', 'AAG', 'AAT', 'AC', 'ACA', 'ACC', 'ACG', 'ACT', 'AG', 'AGA', 
             'AGC', 'AGG', 'AGT', 'AT', 'ATA','ATC', 'ATG', 'ATT', 'CA', 'CAA', 'CAC', 'CAG',
             'CAT', 'CC', 'CCA', 'CCC', 'CCG', 'CCT', 'CG', 'CGA', 'CGC', 'CGG','CGT', 'CT', 
             'CTA', 'CTC', 'CTG', 'CTT', 'GA', 'GAA', 'GAC', 'GAG', 'GAT', 'GC', 'GCA', 'GCC', 
             'GCG', 'GCT', 'GG','GGA', 'GGC', 'GGG', 'GGT', 'GT', 'GTA', 'GTC', 'GTG', 'GTT',
             'TA', 'TAA', 'TAC', 'TAG', 'TAT', 'TC', 'TCA', 'TCC','TCG', 'TCT', 'TG', 'TGA',
             'TGC', 'TGG', 'TGT', 'TT', 'TTA', 'TTC', 'TTG', 'TTT']
    ctd_dict = {x:seq.count(x)*100/n for x in names}
    dis = distance(seq,ctd_dict['A']*n/100.0,ctd_dict['T']*n/100.0,ctd_dict['G']*n/100.0,ctd_dict['C']*n/100.0)
    ctd_dict.update(dis)
    return ctd_dict


def distance(seq,num_A,num_T,num_G,num_C):
    n = len(seq)
    a,t,g,c=0,0,0,0
    A0_dis,A1_dis,A2_dis,A3_dis,A4_dis=0.0,0.0,0.0,0.0,0.0
    T0_dis,T1_dis,T2_dis,T3_dis,T4_dis=0.0,0.0,0.0,0.0,0.0
    G0_dis,G1_dis,G2_dis,G3_dis,G4_dis=0.0,0.0,0.0,0.0,0.0
    C0_dis,C1_dis,C2_dis,C3_dis,C4_dis=0.0,0.0,0.0,0.0,0.0
    for i in range(n):
        if seq[i]=='A':
            a=a+1
            if a == 1:
                A0_dis=(i+1)*100.0/n
            if a == int(round(num_A/4.0)):
                A1_dis=(i+1)*100.0/n
            if a == int(round(num_A/2.0)):
                A2_dis=(i+1)*100.0/n
            if a == int(round((num_A*3/4.0))):
                A3_dis=(i+1)*100.0/n
            if a == num_A:
                A4_dis=(i+1)*100.0/n
        if seq[i]=='T':
            t=t+1
            if t == 1:
                T0_dis=(i+1)*100.0/n
            if t == int(round(num_T/4.0)):
                T1_dis=(i+1)*100.0/n
            if t == int(round((num_T/2.0))):
                T2_dis=(i+1)*100.0/n
            if t == int(round((num_T*3/4.0))):
                T3_dis=(i+1)*100.0/n
            if t == num_T:
                T4_dis=(i+1)*100.0/n
        if seq[i]=='G':
            g=g+1
            if g == 1:
                G0_dis=(i+1)*100.0/n
            if g == int(round(num_G/4.0)):
                G1_dis=(i+1)*100.0/n
            if g == int(round(num_G/2.0)):
                G2_dis=(i+1)*100.0/n
            if g == int(round(num_G*3/4.0)):
                G3_dis=(i+1)*100.0/n
            if g == num_G:
                G4_dis=(i+1)*100.0/n
        if seq[i]=='C':
            c=c+1
            if c == 1:
                C0_dis=(i+1)*100.0/n
            if c == int(round(num_C/4.0)):
                C1_dis=(i+1)*100.0/n
            if c == int(round(num_C/2.0)):
                C2_dis=(i+1)*100.0/n
            if c == int(round(num_C*3/4.0)):
                C3_dis=(i+1)*100.0/n
            if c == num_C:
                C4_dis=(i+1)*100.0/n
    dis = {'A0':A0_dis,'A1':A1_dis,'A2':A2_dis,'A3':A3_dis,'A4':A4_dis,'T0':T0_dis,'T1':T1_dis,'T2':T2_dis,'T3':T3_dis,'T4':T4_dis,'G0':G0_dis,'G1':G1_dis,'G2':G2_dis,'G3':G3_dis,'G4':G4_dis,'C0':C0_dis,'C1':C1_dis,'C2':C2_dis,'C3':C3_dis,'C4':C4_dis}
    return dis

