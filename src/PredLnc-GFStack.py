
import numpy as np
import pandas as pd
import os
import random
import time
from argparse import ArgumentParser
from Bio import SeqIO 
from get_features_module import ProtParam as PP
from get_features_module import fickett
from get_features_module import FrameKmer
from get_features_module import Get_ORF_features as orf
from Bio.SeqUtils import ProtParam
from Bio.SeqUtils import GC
from get_features_module.CTD import CTD
from get_features_module.SNR import SNR
from get_features_module.LncADeepGetFeature import LncADeepGetFeature
from sklearn.ensemble import RandomForestClassifier
from sklearn.externals import joblib
import warnings
warnings.filterwarnings('ignore')

def main():
	parser = ArgumentParser(
			description='PredLnc-GFStack: a global sequence feature-based on Stacked Ensemble learning method for predicting lncRNAs from transcripts.'
						' Liu et al (2019).')
	group = parser.add_argument_group("Method Paramenters")
	group.add_argument('--input', nargs=1, dest='input', type = str, help = "input sequence in fasta format(required)")
	group.add_argument('--model', nargs=1, dest='model', type = str, help = "use [human] model or [mouse] model")
	group.add_argument('--output', nargs=1, dest='output', type = str, help = "output file format like [Sequence Class]")
	args = parser.parse_args()

	if args.input and args.output and args.model:
		input_file = args.input[0]
		output_file = args.output[0]
		model = args.model[0]
		print("[INFO] Getting all features from transcripts...")
		seq_id = get_all_features(input_file,model)
		print("[INFO] Done!")
		print("[INFO] Predicting...")
		label = predict(model)
		print("[INFO] Done!")
		print("[INFO] Saving labels...")
		save_result(label,output_file,seq_id)
		print("[INFO] Finish!")
	else:
		print("Please check your --input and --output and --model parameters!")
	return 0


def GC1(mRNA):
	if len(mRNA) < 3:
		numGC = 0
		mRNA = 'ATG'
	else:
		numGC = mRNA[0::3].count("C") + mRNA[0::3].count("G")
	return numGC*1.0/len(mRNA)*3

def GC2(mRNA):
	if len(mRNA) < 3:
		numGC = 0
		mRNA = 'ATG'
	else:
		numGC = mRNA[1::3].count("C") + mRNA[1::3].count("G")
	return numGC*1.0/len(mRNA)*3

def GC3(mRNA):
	if len(mRNA) < 3:
		numGC = 0
		mRNA = 'ATG'
	else:
		numGC = mRNA[2::3].count("C") + mRNA[2::3].count("G")
	return numGC*1.0/len(mRNA)*3
def Get_ORF_features(input_file):
	transcript_length,first_orf_len,first_orf_Rlen,longest_orf_len,longest_orf_Rlen,\
	integrity,ORF_frame_score,reading_frame_1,reading_frame_2,reading_frame_3,ORF_kmers = [],[],[],[],[],[],[],[],[],[],[]
	seq_id = []
	names = ['A','G','C','T','AA', 'AAA', 'AAC', 'AAG', 'AAT', 'AC', 'ACA', 'ACC', 'ACG', 'ACT', 'AG', 'AGA', 
				 'AGC', 'AGG', 'AGT', 'AT', 'ATA','ATC', 'ATG', 'ATT', 'CA', 'CAA', 'CAC', 'CAG',
				 'CAT', 'CC', 'CCA', 'CCC', 'CCG', 'CCT', 'CG', 'CGA', 'CGC', 'CGG','CGT', 'CT', 
				 'CTA', 'CTC', 'CTG', 'CTT', 'GA', 'GAA', 'GAC', 'GAG', 'GAT', 'GC', 'GCA', 'GCC', 
				 'GCG', 'GCT', 'GG','GGA', 'GGC', 'GGG', 'GGT', 'GT', 'GTA', 'GTC', 'GTG', 'GTT',
				 'TA', 'TAA', 'TAC', 'TAG', 'TAT', 'TC', 'TCA', 'TCC','TCG', 'TCT', 'TG', 'TGA',
				 'TGC', 'TGG', 'TGT', 'TT', 'TTA', 'TTC', 'TTG', 'TTT']
	for seq in SeqIO.parse(input_file,'fasta'):
		a,b,c,d,e,f,g,h,i,j,k = orf.ORF_features(seq.seq)
		transcript_length.append(a);first_orf_len.append(b);first_orf_Rlen.append(c)
		longest_orf_len.append(d);longest_orf_Rlen.append(e);integrity.append(f)
		ORF_frame_score.append(g);reading_frame_1.append(h),reading_frame_2.append(i),reading_frame_3.append(j)
		ORF_kmers.append(k)

	for i in range(len(names)):
		names[i] = str("ORF_kmers_"+names[i])
	ORF_kmers_result = pd.DataFrame(ORF_kmers,columns=names)
	ORF_kmers_result.to_csv("ORF_kmers_features",index=None,sep="\t")
		
	GC1_frame_score,GC2_frame_score,GC3_frame_score,pI_Mw_frame_score=[],[],[],[]
	for i in range(len(reading_frame_1)):
		a = GC1(reading_frame_1[i]);b = GC1(reading_frame_2[i]);c = GC1(reading_frame_3[i])
		GC1_frame_score.append((a - b)**2 + (a - c)**2 +(b- c)**2/2)
		a = GC2(reading_frame_1[i]);b = GC2(reading_frame_2[i]);c = GC2(reading_frame_3[i])
		GC2_frame_score.append((a - b)**2 + (a - c)**2 +(b- c)**2/2)
		a = GC3(reading_frame_1[i]);b = GC3(reading_frame_2[i]);c = GC3(reading_frame_3[i])
		GC3_frame_score.append((a - b)**2 + (a - c)**2 +(b- c)**2/2)
		tmp1 = PP.param(reading_frame_1[i])
		tmp2 = PP.param(reading_frame_2[i])
		tmp3 = PP.param(reading_frame_3[i])	
		pI_Mw_frame_score.append((tmp1[4] - tmp2[4])**2 + (tmp1[4] - tmp3[4])**2 +(tmp2[4]- tmp3[4])/2)

	tmp = {"transcript_length":transcript_length,"first_orf_len":first_orf_len,"first_orf_Rlen":first_orf_Rlen,
	"longest_orf_len":longest_orf_len,"longest_orf_Rlen":longest_orf_Rlen,"integrity":integrity,"ORF_frame_score":ORF_frame_score,
	"GC1_frame_score":GC1_frame_score,"GC2_frame_score":GC2_frame_score,"GC3_frame_score":GC3_frame_score,
	"pI_Mw_frame_score":pI_Mw_frame_score}
	result = pd.DataFrame(tmp)
	result.to_csv("ORF_features",index=None,sep="\t")
	return transcript_length

def coding_nocoding_potential(input_file):
	coding={}
	noncoding={}
	for line in open(input_file).readlines():
		fields = line.split()
		if fields[0] == 'hexamer':continue # pass header
		coding[fields[0]] = float(fields[1])
		noncoding[fields[0]] =  float(fields[2])
	return coding,noncoding
def get_sequence_features(input_file,model):
	if model == "mouse":
		hex_file = "./Mouse_model/Mouse_features_Hexamer.tsv"
	else:
		hex_file = "./Human_model/Human_features_Hexamer.tsv"
	coding,noncoding = coding_nocoding_potential(hex_file)
	hexamer,Instability_index,pI,Gravy,fickett_score,\
	STOP_Codon_Count,STOP_Codon_Frequency,PI_Mw,Mw,\
	transcript_GC1,transcript_GC2,transcript_GC3,transcript_GC,snr,seq_id = [],[],[],[],[],[],[],[],[],[],[],[],[],[],[]
	for seq in SeqIO.parse(input_file,'fasta'):
		transcript_GC.append(GC(seq.seq))
		transcript_GC1.append(GC1(seq.seq))
		transcript_GC2.append(GC2(seq.seq))
		transcript_GC3.append(GC3(seq.seq))
		snr.append(SNR(seq))
		seq_id.append(str(seq.id))
		tmp = PP.param(seq.seq)
		Mw.append(tmp[3])
		Instability_index.append(tmp[0]);pI.append(tmp[1]);Gravy.append(tmp[2])
		PI_Mw.append(tmp[4])
		a = fickett.fickett_value(seq.seq)
		fickett_score.append(a)
		a = FrameKmer.kmer_ratio(seq.seq,6,3,coding,noncoding)
		hexamer.append(a)
		a = (seq.seq).translate().count("*")
		STOP_Codon_Count.append(a)
		STOP_Codon_Frequency.append(float(a)/len(seq.seq))
		
	tmp = {"hexamer_score":hexamer,"Instability_index":Instability_index,
	"pI":pI,"Gravy":Gravy,"fickett_score":fickett_score,"Mw":Mw,
	"PI_Mw":PI_Mw,"STOP_Codon_Count":STOP_Codon_Count,"STOP_Codon_Frequency":STOP_Codon_Frequency,
	"transcript_GC1":transcript_GC1,"transcript_GC2":transcript_GC2,"transcript_GC3":transcript_GC3,"transcript_GC":transcript_GC,"snr":snr}
	result = pd.DataFrame(tmp)
	result.to_csv("sequence_features",index=None,sep="\t")
	return seq_id

def get_special_features(input_file,transcript_length):
	#Get txCdsPredict score, CDS length and CDS percentage
	tmp_cds = "tmp.cds"
	tmp_fa = "tmp.fa"
	CDS_length=[];CDS_percentage=[];txCdsPredict_score=[]

	for seq in SeqIO.parse(input_file,'fasta'):
			 SeqIO.write(seq,tmp_fa,"fasta")
			 os.system("txCdsPredict "+tmp_fa+" "+tmp_cds)
			 if(os.path.getsize(tmp_cds))==0:
					 CDS_length.append(0);txCdsPredict_score.append(0)
			 else:
					 for line in open(tmp_cds).readlines():
							 line = line.split("\t")
							 CDS_length.append(float(line[2])-float(line[1]))
							 txCdsPredict_score.append(float(line[5]))
	for i in range(len(CDS_length)):
			 tmp = float(CDS_length[i]/transcript_length[i])
			 CDS_percentage.append(tmp)
	os.system("rm -f "+tmp_cds)
	os.system("rm -f "+tmp_fa)

	tmp = {"CDS_length":CDS_length,
	"CDS_percentage":CDS_percentage,"txCdsPredict_score":txCdsPredict_score}
	result = pd.DataFrame(tmp)
	result.to_csv("special_features",index=None,sep="\t")

def get_CTD_LncADeep_features(input_file):
	tmp_dict = {}
	ctd = {}
	lnc = {}
	for record in SeqIO.parse(input_file,"fasta"):	
		seq = record.seq
		name = str(record.id).lower()
		tmp_dict[name] = CTD(seq)
		tmp_dict[name].update(LncADeepGetFeature(seq))

	tmp1 = pd.DataFrame(tmp_dict)
	tmp2 = tmp1.T
	tmp2.to_csv("CTD_LncADeep_features.csv",index=None)

#Make table of all results 
def get_all_features(input_file,model):
	transcript_length = Get_ORF_features(input_file)
	seq_id = get_sequence_features(input_file,model)
	get_special_features(input_file,transcript_length)
	get_CTD_LncADeep_features(input_file)
	a = pd.read_csv("ORF_kmers_features",sep="\t")
	b = pd.read_csv("special_features",sep="\t")
	c = pd.read_csv("sequence_features",sep="\t")
	d = pd.read_csv("ORF_features",sep="\t")
	e = pd.read_csv("CTD_LncADeep_features.csv",sep=",")
	os.system("rm -f ORF_kmers_features")
	os.system("rm -f special_features")
	os.system("rm -f sequence_features")
	os.system("rm -f ORF_features")
	os.system("rm -f CTD_LncADeep_features.csv")
	all_X_features = pd.concat([a,b,c,d,e],axis=1)
	all_X_features.to_csv("All_features",index=None,sep="\t")
	return seq_id
	
def predict(model):
	All_proba = []
	label = []
	Selected_features = []
	if model == 'mouse':
			temp = "./Mouse_model/"
	else:
			temp = "./Human_model/"
	filename = temp+model+"_selected_features.csv"
	with open(filename,'r') as f:
		for line in f:
			Selected_features.append(list(line.split("\t"))[:-1])
	for i in range(10):
		All_X = pd.read_csv('All_features',sep="\t")
		All_x = All_X[Selected_features[i]]
		#clf
		if model == 'mouse':
			temp = "./Mouse_model/"
		else:
			temp = "./Human_model/"
		model_name = temp+'RF_model_'+str(i)
		clf = joblib.load(model_name)	
		proba = clf.predict_proba(All_x)[:,1]
		All_proba.append(proba)

	All_proba = np.array(All_proba)
	final_proba = list(np.mean(All_proba, axis=0))
	for i in range(len(final_proba)):
		if final_proba[i] >= 0.5:
			label.append(1)
		else:
			label.append(0)
	os.system("rm -f All_features")
	return label

def save_result(label,output_file,seq_id):
	with open(output_file,'w') as out:
		out.write("Sequence"+'\t'+"Class"+'\n')
		count = 0
		for i in label:
			if i == 1:
				out.write(seq_id[count]+'\t'+'lncRNA'+'\n')
			else:
				out.write(seq_id[count]+'\t'+'pct'+'\n')
			count += 1
			
if __name__ == "__main__":
	main()


