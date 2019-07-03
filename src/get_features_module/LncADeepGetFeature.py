import numpy as np
from .EdpFea import *

#EDP_fea       20
#Kmer_EDP_fea  16
#UTR5_fea      2
#UTR3_fea      2
#fickett_fea   8

def LncADeepGetFeature(seq):
    # search ORF
    seq = seq.strip()
    ORF, UTR5, UTR3 = GetORF_UTR(seq)
    transcript_len = len(seq)

    # The EDP of ORF
    
    if len(ORF) < 6:
        EDP_fea = GetEDP_noORF()
    else:
        EDP_fea = GetEDP(ORF, transcript_len)
   

   
    Kmer_EDP_fea = GetKmerEDP(ORF)

    # UTR feature

    UTR5_fea = str(len(UTR5) * 1.0 / len(seq)) + "\t" + str(GetGC_Content(UTR5))
    UTR3_fea = str(len(UTR3) * 1.0 / len(seq)) + "\t" + str(GetGC_Content(UTR3))
    

    # Fickett feature
    
    A_pos_fea = GetBasePositionScore(seq, 'A')
    C_pos_fea = GetBasePositionScore(seq, 'C')
    G_pos_fea = GetBasePositionScore(seq, 'G')
    T_pos_fea = GetBasePositionScore(seq, 'T')
    base_ratio = GetBaseRatio(seq)

    fickett_fea = "\t".join([str(A_pos_fea), str(C_pos_fea), str(G_pos_fea), str(T_pos_fea)])
    feature = "\t".join([EDP_fea, Kmer_EDP_fea, UTR5_fea, UTR3_fea, fickett_fea])
    feature_name = ['EDP_fea_A', 'EDP_fea_C', 'EDP_fea_D', 'EDP_fea_E', 'EDP_fea_F',
       'EDP_fea_G', 'EDP_fea_H', 'EDP_fea_I', 'EDP_fea_K', 'EDP_fea_L',
       'EDP_fea_M', 'EDP_fea_N', 'EDP_fea_P', 'EDP_fea_Q',
       'EDP_fea_R', 'EDP_fea_S', 'EDP_fea_T', 'EDP_fea_V',
       'EDP_fea_W', 'EDP_fea_Y','Kmer_EDP_ORF_AA', 'Kmer_EDP_ORF_AC', 'Kmer_EDP_ORF_AG',
       'Kmer_EDP_ORF_AT', 'Kmer_EDP_ORF_CA', 'Kmer_EDP_ORF_CC',
       'Kmer_EDP_ORF_CG', 'Kmer_EDP_ORF_CT', 'Kmer_EDP_ORF_GA',
       'Kmer_EDP_ORF_GC', 'Kmer_EDP_ORF_GG', 'Kmer_EDP_ORF_GT',
       'Kmer_EDP_ORF_TA', 'Kmer_EDP_ORF_TC', 'Kmer_EDP_ORF_TG',
       'Kmer_EDP_ORF_TT','UTR5_perc_len', 'UTR5_GC_content','UTR3_perc_len', 'UTR3_GC_content','A_pos_fickett', 
        'C_pos_fickett', 'G_pos_fickett', 'T_pos_fickett']
    feature_dict = {}
    feature_data = np.fromstring(feature,sep=' ')
    if len(feature_data)==len(feature_name):
        for i in range(len(feature_data)):
            feature_dict[feature_name[i]]=feature_data[i]
    return feature_dict
    

