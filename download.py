from Bio import Entrez
from Bio import SeqIO
import os,argparse
parser = argparse.ArgumentParser(description='This script was used to download genbank file')
parser.add_argument('-i','--id',help='Please input sgRNA file,file requires two columns, the first column is id, the second column is sequence, and no parentheses can exist.',required=True)
parser.add_argument('-o','--out_dir',help='Please input	 out_put directory path;default cwd',default = os.getcwd(),required=False)
args = parser.parse_args()
if os.path.exists(args.out_dir):
	kout=os.path.abspath(args.out_dir)
else:
	os.mkdir(args.out_dir)
	kout=os.path.abspath(args.out_dir)
def detective(path):
	has_name=[]
	namelst=os.listdir(path)
	for name in namelst:
		has_name.append(name.split(".")[0])
	return has_name
def genbank_download(genid):
	Entrez.email = "liangxue@genomics.cn"
	handle = Entrez.esummary(db='gene',id=genid)
	try:
		records = Entrez.read(handle)
	except:
		print("this id can not analyse")
	else:
		name=records['DocumentSummarySet']['DocumentSummary'][0]['Name']
		ids=records['DocumentSummarySet']['DocumentSummary'][0]['GenomicInfo'][0]['ChrAccVer']
		start=records['DocumentSummarySet']['DocumentSummary'][0]['GenomicInfo'][0]['ChrStart']
		stop=records['DocumentSummarySet']['DocumentSummary'][0]['GenomicInfo'][0]['ChrStop']
		if int(start)<int(stop):
			hd = Entrez.efetch(db="nucleotide",id=ids,rettype='gb',retmode='text',seq_start=start,seq_stop=stop)
		else:
			hd = Entrez.efetch(db="nucleotide",id=ids,rettype='gb',retmode='text',seq_start=stop,seq_stop=start)
		seq = SeqIO.read(hd,'gb')
		path=kout+"/%s.gb"%name
		fw = open(path,'w')
		SeqIO.write(seq,fw,'gb')
		fw.close()
		os.getcwd()
	return genid


if __name__ == '__main__':
	path1=os.path.abspath(args.id)
	file1=open(path1)
	id_lst=file1.readlines()
	id_lst=[genid.strip("\n") for genid in id_lst]
	idlst=[]
	try:
		with open('log.txt') as k:
			for line in k:
				idlst.append(line.strip())
	except:
		print("new down")
	for genid in id_lst:
		if genid in idlst:
			print("alread exits")
			continue
		else:
			with open('log.txt','a') as f:
				f.write(genid+"\n")
			genbank_download(genid)
	print("gb files have been download ")
