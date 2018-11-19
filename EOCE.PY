from Bio import Entrez
from Bio import SeqIO
import os,argparse,sys,glob,os.path,re
parser = argparse.ArgumentParser(description='This script was used to download genbank file')
parser.add_argument('-o','--out_dir',help='Please input	 out_put directory path;default cwd',default = os.getcwd(),required=False)
parser.add_argument('-u','--out_dir2',help='Please input  out_put directory path;default cwd',default = os.getcwd(),required=False)
parser.add_argument('-s','--out_dir3',help='Please input  out_put directory path;default cwd',default = os.getcwd(),required=False)
args = parser.parse_args()
#批量下载genbank序列^^^
def detective(path):
	has_name=[]
	namelst=os.listdir(path)
	for name in namelst:
		has_name.append(name.split(".")[0])
	return has_name
input_name=[]
if os.path.exists(args.out_dir):
	kout=os.path.abspath(args.out_dir)
else:
	os.mkdir(args.out_dir)
	kout=os.path.abspath(args.out_dir)
if os.path.exists(args.out_dir2):
	dout=os.path.abspath(args.out_dir2)
else:
	os.mkdir(args.out_dir2)
	dout=os.path.abspath(args.out_dir2)
if os.path.exists(args.out_dir3):
	sout=os.path.abspath(args.out_dir3)
else:
	os.mkdir(args.out_dir3)
	sout=os.path.abspath(args.out_dir3)
gene = {}
input = os.listdir(kout)
for line in input :
	line = line.strip(".gb")
	if line in detective(sout):
		print("already exits")
		continue
	else:
		input_name.append(line)
		gene[line] = line
genbank=[genbank+".gb" for genbank in input_name]
quan_genbank=[os.path.join(kout,haha) for haha in genbank]
for gdkfile in quan_genbank :
	name = os.path.basename(gdkfile)
	input_handle  = open(gdkfile, "r")
	gene_file = dout+'/'+name+".genome.fa"
	gene_handle = open(gene_file, "w")
	cds_locat_file = dout+'/'+name+".cds_location.txt"
	cds_locat_handle = open(cds_locat_file, "w")
	for seq_record in SeqIO.parse(input_handle, "genbank") :
		print("Dealing with GenBank record %s" % seq_record.id)
		gene_handle.write(">%s %s\n%s\n" % (seq_record.id,seq_record.description,seq_record.seq))
		for seq_feature in seq_record.features :
			geneSeq = seq_feature.extract(seq_record.seq)
			if seq_feature.type=="CDS" :
				#assert len(seq_feature.qualifiers['translation'])==1
				if seq_feature.qualifiers['gene'][0] in gene.keys():
					cds_locat_handle.write(">%s location %s\n" % (seq_feature.qualifiers['gene'][0],seq_feature.location))
					#print(seq_feature.qualifiers['translation'],seq_feature.location)
gene_handle.close()
cds_locat_handle.close()
input_handle.close()
#genePEP.close()
#geneCDS.close()
print("cdss had been catch")
#提取出cds的位置^^^
#该函数将提取所有的位置信息，并返回【‘起始位置，结束位置’……】这样的一个列表
def lst_tiqu(str_seq):
	r=re.compile('\d+:\d+')
	lst=r.findall(str_seq)
	lst=[lst.replace(":",",") for lst in lst]
	return lst
#该函数将上一步形成的列表改变为【【’起始位置’，‘结束位置’】……】这样的嵌套列表
def pinjie(lst):
	lst=lst.split(",")
	return lst
#该函数将两个列表合并，并取列表的中间两位数并重新组成一个包含起始位置和结束位置的字符串
def internal(lst_A,lst_B):
	lst_A.extend(lst_B)
	lst=[int(x) for x in lst_A]
	lst.sort()
	new_lst="%d,%d"%(lst[1],lst[2])
	return new_lst
#该函数提取位置重叠区域，返回重叠区域列表
def jiaoji(lst_A,lst_B):
	new_lst=[]
	for lst1 in lst_A:
		for lst2 in lst_B:
			if int(lst2.split(",")[0])< int(lst1.split(",")[1]) and int(lst2.split(",")[1])> int(lst1.split(",")[0]):
				new_lst.append(internal(pinjie(lst1),pinjie(lst2)))
			else:
				pass
	return new_lst
#将所有的基因名字组成字符串为input_name
with open("./positions.txt",'w') as contig_file:
	for gene_name in input_name:
		new_lst=[]
		cdsloc_file=""
		locpath=""
		locpath=dout+"/"+gene_name+".gb.cds_location.txt"
		cdsloc_file=open(locpath)
		str_lst=cdsloc_file.readlines()
		if len(str_lst)==1:
			new_lst=lst_tiqu(str_lst[0])
			print("we are caculating %s sequences for 1"%gene_name)
			if len(new_lst)<=2:
				contig_file.write(gene_name+"\t"+str(new_lst)+"\n")
			else:
				contig_file.write(gene_name+"\t"+str(new_lst[1:-1])+"\n")
		elif len(str_lst)==2:
			new_lst=jiaoji(lst_tiqu(str_lst[0]),lst_tiqu(str_lst[1]))#先第一个和第二个相互求交集
			print("we are caculating %s sequences  for 2"%gene_name)
			if len(new_lst)<=2:
				contig_file.write(gene_name+"\t"+str(new_lst)+"\n")
			else:
				contig_file.write(gene_name+"\t"+str(new_lst[1:-1])+"\n")
		else:
			for i in range(len(str_lst)):
				if i==0:
					new_lst=jiaoji(lst_tiqu(str_lst[i]),lst_tiqu(str_lst[i+1]))
				elif i==1:
					pass
				else:
					new_lst=jiaoji(new_lst,lst_tiqu(str_lst[i]))#从第三个开始，依次计算前一个结果和后一个位置的交集
			print("we are caculating %s sequences for 3"%gene_name)
			if len(new_lst)<=2:
				contig_file.write(gene_name+"\t"+str(new_lst)+"\n")
			else:
				contig_file.write(gene_name+"\t"+str(new_lst[1:-1])+"\n")
	print("done")
#以上函数提取所有cds序列的交集并输出到文件
with open("./positions.txt") as f:
	str1_lst= f.readlines()
#该函数将位置信息解析为[[start,stop],....]
def pos_parse(str1):
	hah=[]
	r=re.compile('\d+,\d+')
	lst=r.findall(str1)
	new_lst=[lst.split(",") for lst in lst]
	for lst in new_lst:
		hah.append([int(nn) for nn in lst])
	return hah
#该函数打开以位置信息命名的基因名字为名字的fasta序列
def opfile(str1):
	gene_name=str1.split("\t")[0]#提取位置信息中基因的名字
	genome_path=dout+"/%s.gb.genome.fa"%gene_name
	geno_file=open(genome_path)
	seq= geno_file.readlines()[1]
	geno_file.close()
	return seq
#提取位置信息中基因的名字
def get_name(str1):
	gene_name=str1.split("\t")[0]
	return gene_name
seq_lst=[]
for str1 in str1_lst:
	seq_lst.append([opfile(str1)[hahh[0]:hahh[1]] for hahh in pos_parse(str1)])

for i in range(len(seq_lst)):##嵌套循环将所有基因的所有序列依次写入文件
	print(">","%s"%get_name(str1_lst[i]))
	print("\n")
	seq_path=sout+"/%s.cds_coverage.fa"%get_name(str1_lst[i])
	s=open(seq_path,'w')#写入文件
	for j in range(len(seq_lst[i])):
		print(seq_lst[i][j])
		print("\n")
		s.write(">"+"%s"%get_name(str1_lst[i])+"(%s--%s)"%(pos_parse(str1_lst[i])[j][0],pos_parse(str1_lst[i])[j][1])+"\n")
		s.write(seq_lst[i][j]+"\n")
s.close()









