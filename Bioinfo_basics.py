#1. find multiple genes 
"""
"""
genome_file = 'C:/Users/alejc/Desktop/' #main file 
genelistfile = 'C:/Users/alejc/Desktop/' #genes looking to be found

from Bio import SeqIO
with open (genelistfile,'r') as input_file:
    gene_names=[line.strip('\n') for line in input_file]

len(gene_names)
gene_names[0:5] #obtain gene names

gb_object=SeqIO.read(genome_file,'gb')
allgenes=[feature for feature in gb_object.features if feature.type =='gene']
len(allgenes)

gene_sequences=[]
for gene in allgenes: #check if genes are present
    if 'gene' in gene.qualifiers.keys():
        gene_name=gene.qualifiers['gene'][0]
        if gene_name in gene_names:
            extract=gene.extract(gb_object)
            extract.id=gene_name
            extract.description=''
            gene_sequences.append(extract)
            print('gene %s has been found'%gene_name)

len(gene_sequences)

#2. DNA complementer 
"""
"""
#complementary ADN sequences code 
dna = 'ATCTAGGATATAC'
#create dictionary
base_complements={'A':'T','T':'A','C':'G','G':'C'}
#complements of seq
s_c=''
for base in dna:
    base_complement=base_complements[base]
    s_c+=base_complement
  
print(s_c) #checking

#functions + loops
dnaseq = ['ATCGA', 'TTAGC', 'ACCTAG']
s_c=''
def sequence_complement_finder(seq):
    base_complements={'A':'T','T':'A','C':'G','G':'C'}
    s_c= [base_complements[base] for base in seq]
    s_c= ''.join(s_c)
    return s_c

for seq in dnaseq:
    s_c=sequence_complement_finder(seq)
    print('original sequence:', seq)
    print('sequence complement:', s_c)

complement_sequences=[]
for seq in dnaseq:
    s_c= sequence_complement_finder(seq)
    complement_sequences.append(s_c)
print('original:     ',dnaseq)
print('complement:     ',complement_sequences)

#3. some others-- seq 
a = 992
b = 896
c = a^2 + b^2

#how to find length of sequence
dna1 = 'AGTCCTAGAAGATCCAA'
len(dna1)

#presence or absence of characters
dna = "ATCGGA"
'A' in dna

#number of nucleotides
dna2 = "AGTAGGATCGATTATATTAA"
dna2.count("A")
dna2.count("G")

#percentage
gc = dna2.count("G") + dna2.count("C")
length = len(dna2)
gc_percent = gc/length * 100
print(gc_percent)
           
#combine strigns/sequences
dna1 = 'ATGA'
dna2 = 'GCCG'
com = dna1 + dna2
           
print(com)

#4. Download fasta and genbank files
#nucleotide sequences from ncbi
from Bio import Entrez
from Bio import SeqIO
Entrez.email=''
handle=Entrez.efetch(db='nuccore',id='34577062')
print(handle.read())
#fasta way
handle=Entrez.efetch(db='nuccore',id='34577062',rettype='fasta')
print(handle.read())
#resumed info way
handle=Entrez.efetch(db='nucleotide',id='NM_001126.2',rettype='gb',retmode='text')
print(handle.read())
#analyze fasta
handle=Entrez.efetch(db='nucleotide',id='NM_001126.2',rettype='fasta',retmode='text')
record=SeqIO.read(handle,'fasta')
record.id
record.name
record.description
seq=record.seq
print(seq[0:10])
len(seq)
print(seq)

#ETC
handle=Entrez.efetch(db='nucleotide',id='NM_001126.2',rettype='gb',retmode='text')
record=SeqIO.read(handle,'gb')
record.id
features=record.features
len(features)
print(features)
seq=record.seq
print(len(seq))

#Multiple genomes
handle=Entrez.efetch(db='nuccore',id='34577062 , 186972394',rettype='fasta')
records=SeqIO.parse(handle, 'fasta')
records=[i for i in records]
len(records)
first_record=records[0]
second_record=records[1]
first_record.id
second_record.id

#download and save (gb)
handle=Entrez.efetch(db='nucleotide',id='34577062',rettype='gb')
record=SeqIO.read(handle,'gb')
outputname='C:/Users/alejc/Desktop/kuh/ADSS.gb'
SeqIO.write(record,outputname,'gb')

#fasta
handle=Entrez.efetch(db='nucleotide',id='34577062',rettype='fasta')
record=SeqIO.read(handle,'fasta')
outputname='C:/Users/alejc/Desktop/kuh/ADSS.fasta'
SeqIO.write(record,outputname,'fasta')

#for multiple just change read for parse

#5. Feature counts GB
from Bio import SeqIO
import pandas as pd
from collections import Counter
file_path = ''

genbank_obj=SeqIO.read(file_path,'gb')
all_f_types=[feature.type for feature in genbank_obj.features]
len(all_f_types)

feature_types=set(all_f_types)
print(feature_types)

feature_counts=Counter(all_f_types)
feature_counts.keys()
for key,value in feature_counts.items():
    print(key,value)
del feature_counts['source']
del feature_counts['regulatory']
dataframe=pd.DataFrame(feature_counts.items(),columns=['Feature','Count'])
dataframe
outputfile= 'C:/.csv'
dataframe.to_csv(outputfile,index=False) 

#Feature count MULTIPLE
import glob
from Bio import SeqIO
from collections import Counter
import pandas as pd
import os
directory = 'C:/Users/alejc/Desktop/callitkarma'
gfiles = glob.glob('%s/*.gb'%directory)
gfiles
#function to read and get features genbank
def read_file(gfile):
    genbank_obj=SeqIO.read(gfile,'gb')
    features=genbank_obj.features
    feature_types=[feature.type for feature in features]
    return feature_types

def count_features(feature_types):
    feature_count = Counter(feature_types)
    print('features have been counted')
    return feature_count

#function to identify features in gb files
def scan_all_features(files):
    allfeatures=[]
    for gfile in gfiles:
        feature_types=read_file(gfile)
        allfeatures.extend(feature_types)

        allfeatures = set(allfeatures)
        allfeatures = list(allfeatures)
        print('all features have been identified')
        return allfeatures

allfeatures=scan_all_features(gfiles)
print(allfeatures)
#emptylist just bc
allfeature_count=[]
for gfile in gfiles:
    directory,filename=os.path.split(gfile)
    filename=filename.strip('.gb')
    feature_types=read_file(gfile)
    feature_count=count_features(feature_types)
    temp_count=[]
    temp_count.append(filename)
    for feature in allfeatures:
        if feature in feature_count.keys():
            temp_count.append(feature_count[feature])
        else:
            temp_count.append(0)
    allfeature_count.append(temp_count)

len(allfeature_count)
allfeature_count

#Database
columns=[]
columns.append('File')
columns.extend(allfeatures)
columns
dataframe=pd.DataFrame(allfeature_count,columns = columns)
print(dataframe)
finaldata=dataframe.set_index('File')
finaldata
columns_to_delete=['source']
finaldata.drop(columns=columns_to_delete,inplace=True)
print(finaldata)

#specific location
finaldata.loc['V521',:]
finaldata.loc['V521','gene']

#output file
outputfile= ''
finaldata.to_csv(outputfile,index=True)

#6. Venn Diagrams
import matplotlib.pyplot as plt
from matplotlib_venn import venn2, venn2_unweighted, venn2_circles
items=[80,43,6]
labels=['Set1','Set2']
#diagram
venn2(subsets=items,set_labels=labels,set_colors=('red','green'),alpha=0.7)
plt.savefig('')
plt.title('2 circle venn diagram')
plt.show()
#emptycircles
venn2_circles(subsets=items)
plt.show()
#3 circles
from matplotlib_venn import venn3, venn3_unweighted
items2=[88,43,6,62,16,3,18]
labels2=['Set1','Set2','Set3']
venn3(subsets=items2,set_labels=labels2,alpha=0.4)
plt.show()
venn3_unweighted(subsets=items2,set_labels=labels2,alpha=0.4)
plt.show()