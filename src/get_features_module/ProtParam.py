import sys
import re
from Bio.Seq import Seq
from .ORF import ExtractORF
from Bio.SeqUtils import ProtParam
import numpy as np

def mRNA_translate(mRNA):
	return Seq(mRNA).translate()

def protein_param(putative_seqprot):
	return (putative_seqprot.instability_index(),putative_seqprot.isoelectric_point(),putative_seqprot.gravy(),putative_seqprot.molecular_weight())

def param(seq):
	strinfoAmbiguous = re.compile("X|B|Z|J|U",re.I)
	ptU = re.compile("U",re.I)
	seqRNA = ptU.sub("T",str(seq).strip())
	seqRNA = seqRNA.upper()
	CDS_size1,CDS_integrity,seqCDS= ExtractORF(seqRNA).longest_ORF(start=['ATG'],stop=['TAA','TAG','TGA'])
	seqprot = mRNA_translate(seqCDS)
	pep_len = len(seqprot.strip("*"))
	newseqprot = strinfoAmbiguous.sub("",str(seqprot))
	protparam_obj = ProtParam.ProteinAnalysis(str(newseqprot.strip("*")))
	if pep_len > 0:
		Instability_index,PI,Gravy,Mw = protein_param(protparam_obj)
		pI_Mw = np.log10((float(Mw)/PI) + 1)
	else:
		Instability_index = 0.0
		PI=0.0
		Gravy=0.0
		Mw = 0.0
		pI_Mw = 0.0
	return(Instability_index,PI,Gravy,Mw,pI_Mw)
