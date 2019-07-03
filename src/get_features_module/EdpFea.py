import os,sys
import numpy as np

## AA list 
_AA_list = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y']

## 2mer list
_DNA = ['A', 'C', 'G', 'T']
_Kmer_list = []
for dna1 in _DNA:
    for dna2 in _DNA:
        _Kmer_list.append(dna1+dna2)

## 3mer list
_3mer_list = []
for dna1 in _DNA:
    for dna2 in _DNA:
        for dna3 in _DNA:
            _3mer_list.append(dna1+dna2+dna3)

## IUPAC code
_IUPAC = {'A':'A','C':'C', 'G':'G', 'T':'T', 'R':'AG', 'Y':'CT', 'M':'AC', 'K':'GT', 'S':'CG', 'W':'AT', 'H':'ACT', \
          'B':'CGT', 'V':'ACG', 'D':'AGT', 'N':'ACGT'}

def IUPAC_2mer(seq):
    '''Return a list of all possible 2mers of the sequence'''

    kmer_list = []
    for dna1 in _IUPAC[seq[0]]:
        for dna2 in _IUPAC[seq[1]]:
            kmer_list.append(dna1+dna2)
    return kmer_list

def IUPAC_3mer(seq):
    '''Return a list of all possible 3mers of the sequence'''

    kmer_list = []
    for dna1 in _IUPAC[seq[0]]:
        for dna2 in _IUPAC[seq[1]]:
            for dna3 in _IUPAC[seq[2]]:
                if Codon2AA2(dna1+dna2+dna3) != "J":
                    kmer_list.append(dna1+dna2+dna3)
    return kmer_list

def SixMer2AA(seq):
    '''Convert 6mer to 2 AA'''

    return Codon2AA2( seq[0:3] ) + Codon2AA2( seq[3:6] )

def GetBasePositionScore(seq, base):
    '''compute base(A, C, G, T) position score based
       on Fickett's methods
    '''
    base_num = [0, 0, 0]
    tmp = len(seq) // 3

    for i in range(0, tmp):
        for j in range(0, 3):
            if seq[j + 3*i] == base:
                base_num[j] += 1

    base_pos_score = max(base_num) * 1.0 / (min(base_num) + 1)

    return base_pos_score


def GetBaseRatio(seq):
    '''calculate the A, C, G, T percentage of ORF'''
    A, C, G, T = 1e-9, 1e-9, 1e-9, 1e-9 
    for i in range(0, len(seq)):
        if seq[i] == 'A':
            A += 1
        elif seq[i] == 'C':
            C += 1
        elif seq[i] == 'G':
            G += 1
        elif seq[i] == 'T':
            T += 1

    A = A * 1.0 / (len(seq) + 4e-9)
    C = C * 1.0 / (len(seq) + 4e-9)
    G = G * 1.0 / (len(seq) + 4e-9)
    T = T * 1.0 / (len(seq) + 4e-9)

    return [A, C, G, T]


def GetGC_Content(seq):
    '''calculate GC content of sequence'''
    A, C, G, T = 1e-9, 1e-9, 1e-9, 1e-9 
    for i in range(0, len(seq)):
        if seq[i] == 'A':
            A += 1
        elif seq[i] == 'C':
            C += 1
        elif seq[i] == 'G':
            G += 1
        elif seq[i] == 'T':
            T += 1

    GC = (G + C) / (A + C + G + T)

    return GC 


def GetORF_UTR(seq):
    '''Get ORF and UTR from sequence'''
    STP = {}
    STP[0] = []
    STP[1] = []
    STP[2] = []
    STP[0].append(0)
    STP[1].append(1)
    STP[2].append(2)
    AAnum = len(seq) // 3

    for i in range(0, 3):
        for j in range(0, AAnum):
            tmp = seq[(i+3*j):(i+3*(j+1))]
            if tmp == 'TAG' or tmp == 'TAA' or tmp == 'TGA':
                STP[i].append(i+3*j)

    ORF = {}

    for i in range(0,3):
        if len(STP[i]) < 2:
            continue
        for j in range(1, len(STP[i])):
            tmpN = (STP[i][j] - STP[i][j-1])//3
            for k in range(0, tmpN):
                tmpS = seq[ (STP[i][j-1] + 3*k):(STP[i][j-1] + 3*(k+1)) ]
                if tmpS == 'ATG':
                    ORF[3*k + STP[i][j-1]] = STP[i][j]
                    break

    # longest ORF
    ORFseq = []
    ORFlen = []
    ORFstart = []
    ORFend = []
    for (k,v) in ORF.items():
        ORFseq.append(seq[k:(v-3)])
        ORFlen.append(v - k)
        ORFstart.append(k)
        ORFend.append(v)

    # ORF doesn't exist
    if ORF:
        ORF_l = ORFseq[np.argmax(ORFlen)]
        UTR5 = ''
        UTR3 = ''
        if len(seq[ 0:ORFstart[np.argmax(ORFlen)] ]) > 0:
            UTR5 = seq[ 0:ORFstart[np.argmax(ORFlen)] ]
        if len(seq[ORFend[np.argmax(ORFlen)]:]) > 0:
            UTR3 = seq[ORFend[np.argmax(ORFlen)]:]
        return [ORF_l, UTR5, UTR3]
    else:
        return ['', '', '']



def GetEDP_noORF():
    '''return the default values'''

    Codon = {}
    for aa in _AA_list:
        Codon[aa] = 1e-9 

    sum_codon = 1e-9 * 20 

    H = 0.0
    for (k,v) in Codon.items():
        Codon[k] /= sum_codon
        Codon[k] = -Codon[k] * np.log2(Codon[k])
        H += Codon[k]

    EDP = {}
    for (k,v) in Codon.items():
        EDP[k] = Codon[k] / H

    outline = ''
    for i in range(20):
        outline += str(0) + "\t"

    return outline 


def GetEDP(seq, transcript_len):
    '''get features including: ORF length, ORF ratio, ORF EDP of codon'''

    Codon = {}
    for aa in _AA_list:
        Codon[aa] = 1e-9 

    sum_codon = 1e-9 * 20 

    if(len(seq) > 3):
        num = len(seq) // 3
        for i in range(0,num) :
            if Codon2AA2( seq[i*3:(i+1)*3] ) == "J":
                continue
            ## IUPAC codon
            elif Codon2AA2( seq[i*3:(i+1)*3] ) == "Z":
                tmp_kmer_list = IUPAC_3mer(seq[i*3:(i+1)*3])
                for tmp_kmer in tmp_kmer_list:
                    Codon[ Codon2AA2(tmp_kmer) ] += 1.0 / len(tmp_kmer_list)
                sum_codon += 1.0
            else:
                Codon[ Codon2AA2( seq[i*3:(i+1)*3] ) ] += 1.0 
                sum_codon += 1.0

        H = 0.0
        for (k,v) in Codon.items():
            Codon[k] /= sum_codon
            Codon[k] = -Codon[k] * np.log2(Codon[k])
            H += Codon[k]

        EDP = {}
        for (k,v) in Codon.items():
            EDP[k] = Codon[k] / H

        outline = ''
        for (k,v) in EDP.items():
            outline += str(v) + "\t"
        
        return outline 


def GetKmerEDP(seq):

    Kmer = {}
    for aa in _Kmer_list:
        Kmer[aa] = 1e-9

    sum_Kmer = 1e-9 * 16 

    if(len(seq) > 3):
        for i in range(0,len(seq)-1) :
            ## IUPAC kmer
            if seq[ i:(i+2) ] not in _Kmer_list:
                tmp_kmer_list = IUPAC_2mer(seq[i:(i+2)])
                for tmp_kmer in tmp_kmer_list:
                    Kmer[ tmp_kmer ] += 1.0 / len(tmp_kmer_list)
            else:
                Kmer[ seq[ i:(i+2)] ] += 1.0 
            sum_Kmer += 1.0

        H = 0.0
        for (k,v) in Kmer.items():
            Kmer[k] /= sum_Kmer
            Kmer[k] = -Kmer[k] * np.log2(Kmer[k])
            H += Kmer[k]

        EDP = {}
        for (k,v) in Kmer.items():
            EDP[k] = Kmer[k] / H

        outline = ''
        for (k,v) in EDP.items():
            outline += str(v) + "\t"

        return outline

    else:
        return GetKmerEDP_Default()


def GetKmerEDP_Default():
    '''kmer entropy density'''
    Kmer = {}
    for aa in _Kmer_list:
        Kmer[aa] = 1e-9

    sum_Kmer = 1e-9 * 16 

    H = 0.0
    for (k,v) in Kmer.items():
        Kmer[k] /= sum_Kmer
        Kmer[k] = -Kmer[k] * np.log2(Kmer[k])
        H += Kmer[k]

    EDP = {}
    for (k,v) in Kmer.items():
        EDP[k] = Kmer[k] / H

    outline = ''
    for (k,v) in EDP.items():
        outline += str(v) + "\t"

    return outline

def Codon2AA2(codon):
    '''convert codon to aa'''
    if codon == "TTT" or codon == "TTC":
        return 'F'
    elif codon == 'TTA' or codon == 'TTG' or codon == 'CTT' or codon == 'CTA' or codon == 'CTC' or codon == 'CTG':
        return 'L'
    elif codon == 'ATT' or codon == 'ATC' or codon == 'ATA':
        return 'I'
    elif codon == 'ATG':
        return 'M'
    elif codon == 'GTA' or codon == 'GTC' or codon == 'GTG' or codon == 'GTT':
        return 'V'
    elif codon == 'GAT' or codon == 'GAC':
        return 'D'
    elif codon == 'GAA' or codon == 'GAG':
        return 'E'
    elif codon == 'TCA' or codon == 'TCC' or codon == 'TCG' or codon == 'TCT':
        return 'S'
    elif codon == 'CCA' or codon == 'CCC' or codon == 'CCG' or codon == 'CCT':
        return 'P'
    elif codon == 'ACA' or codon == 'ACG' or codon == 'ACT' or codon == 'ACC':
        return 'T'
    elif codon == 'GCA' or codon == 'GCC' or codon == 'GCG' or codon == 'GCT':
        return 'A'
    elif codon == 'TAT' or codon == 'TAC':
        return 'Y'
    elif codon == 'CAT' or codon == 'CAC':
        return 'H'
    elif codon == 'CAA' or codon == 'CAG':
        return 'Q'
    elif codon == 'AAT' or codon == 'AAC':
        return 'N'
    elif codon == 'AAA' or codon == 'AAG':
        return 'K'
    elif codon == 'TGT' or codon == 'TGC':
        return 'C'
    elif codon == 'TGG':
        return 'W'
    elif codon == 'CGA' or codon == 'CGC' or codon == 'CGG' or codon == 'CGT':
        return 'R'
    elif codon == 'AGT' or codon == 'AGC':
        return 'S'
    elif codon == 'AGA' or codon == 'AGG':
        return 'R'
    elif codon == 'GGA' or codon == 'GGC' or codon == 'GGG' or codon == 'GGT':
        return 'G'
    # stop codon
    elif codon == 'TAA' or codon == 'TAG' or codon == 'TGA':
        return 'J'
    else:
        return 'Z'     ## IUPAC Ambiguity Codes 


