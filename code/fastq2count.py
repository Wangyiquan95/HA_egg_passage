#!/usr/bin/python
import sys
from Bio import SeqIO
from collections import Counter
import os

def translation(seq):
  dnamap = {"TTT":"F", "TTC":"F", "TTA":"L", "TTG":"L",
    "TCT":"S", "TCC":"S", "TCA":"S", "TCG":"S",
    "TAT":"Y", "TAC":"Y", "TAA":"_", "TAG":"_",
    "TGT":"C", "TGC":"C", "TGA":"_", "TGG":"W",
    "CTT":"L", "CTC":"L", "CTA":"L", "CTG":"L",
    "CCT":"P", "CCC":"P", "CCA":"P", "CCG":"P",
    "CAT":"H", "CAC":"H", "CAA":"Q", "CAG":"Q",
    "CGT":"R", "CGC":"R", "CGA":"R", "CGG":"R",
    "ATT":"I", "ATC":"I", "ATA":"I", "ATG":"M",
    "ACT":"T", "ACC":"T", "ACA":"T", "ACG":"T",
    "AAT":"N", "AAC":"N", "AAA":"K", "AAG":"K",
    "AGT":"S", "AGC":"S", "AGA":"R", "AGG":"R",
    "GTT":"V", "GTC":"V", "GTA":"V", "GTG":"V",
    "GCT":"A", "GCC":"A", "GCA":"A", "GCG":"A",
    "GAT":"D", "GAC":"D", "GAA":"E", "GAG":"E",
    "GGT":"G", "GGC":"G", "GGA":"G", "GGG":"G",
    "---":"-"}
  pep = []
  i = 0
  while i < len(seq):
    codon = seq[i:i + 3].upper()
    aa = dnamap[codon]
    pep.append(aa)
    i = i + 3
  pep = ''.join(pep)
  return pep

def sum_mut(aa1,aa2):
    return sum ( aa1[i] != aa2[i] for i in range(len(aa1)) )


def call_mutid(mutpep,refseq,shift):
  mut_id_ls = []
  assert (len(mutpep) == len(refseq))
  for n in range(len(mutpep)):
    pos = n+shift
    if refseq[n]!=mutpep[n]:
       mut_id_ls.append(refseq[n]+str(pos)+mutpep[n])
  return mut_id_ls

def Call_barcode_type(barcode,fqfilename):
    class1 = ["KSGP-SGGP-01", "KSGP-SGGP-02", "KSGP-SGGP-03", "KSGP-SGGP-04", "KSGP-SGGP-05", "KSGP-SGGP-06",
              "KSGP-SGGP-07", "KSGP-SGGP-08", "KSGP-SGGP-09", "KSGP-SGGP-10", "KSGP-SGGP-11", "KSGP-SGGP-12",
              "KSGP-SGGP-13", "KSGP-SGGP-14", "KSGP-SGGP-15"]
    class2 = ["KSVL-SWGP-SWVL-07", "KSVL-SWGP-SWVL-08", "KSVL-SWGP-SWVL-09"]
    class3 = ["KSVL-SWGP-SGVL-01","KSVL-SWGP-SGVL-02","KSVL-SWGP-SGVL-03",
              "KSVL-SWGP-SGVL-04","KSVL-SWGP-SGVL-05","KSVL-SWGP-SGVL-06"]
    class4 = ["KSVL-SWGP-10","KSVL-SWGP-11","KSVL-SWGP-12","KSVL-SWGP-13","KSVL-SWGP-14","KSVL-SWGP-15"]
    filename = (os.path.basename(fqfilename)).rsplit('.')[0]

    if filename in class1:
        if barcode == 'GTC':
            type = 'KSGP'
        elif barcode == 'CGA':
            type = 'SGGP'
        else:
            type = 'other'
    elif filename in class2:
        if barcode == 'GTC':
            type = 'KSVL'
        elif barcode == 'CGA':
            type = 'SWGP'
        elif barcode == 'TCG':
            type = 'SWVL'
        else:
            type = 'other'
    elif filename in class3:
        if barcode == 'GTC':
            type = 'KSVL'
        elif barcode == 'CGA':
            type = 'SWGP'
        elif barcode == 'TCG':
            type = 'SGVL'
        else:
            type = 'other'
    elif filename in class4:
        if barcode == 'GTC':
            type = 'KSVL'
        elif barcode == 'CGA':
            type = 'SWGP'
        else:
            type = 'other'
    else:
        type = 'other'
    return type
def cal_ref_dict(ref):
    Rrecords = SeqIO.parse(ref, "fasta")
    ref_dict={}
    for record in Rrecords:
        seq = str(record.seq)
        id = str(record.id)
        ref_dict[id]=seq
    return ref_dict

def cal_fastq_dic(fastq,ref_dict):

    print ("reading %s" % fastq)
    Rrecords = SeqIO.parse(fastq,"fastq")
    mut_id_ls = []
    error_read=0
    len_error=0
    shift=110
    for record in Rrecords:
        seq = str(record.seq)
        #filter out ambiguous read
        if 'N' in seq:continue
        barcode = seq[:3]
        seq=seq[5:-2]
        type = Call_barcode_type(barcode,fastq)
        if type == 'other':continue
        ref = ref_dict[type]
        if len(seq) != len(ref):
            len_error+=1
        if len(seq) == len(ref):
            mut_aa = translation(seq)
            ref_aa = translation(ref)
            if sum_mut(mut_aa, ref_aa) == 0:
                mut_id = type+'|'+'WT'
                mut_id_ls.append(mut_id)
            else:
                call_mut_id = call_mutid(mut_aa, ref_aa, shift)
                mut_ids = "-".join(sorted(call_mut_id, key=lambda x: int(x[1:-1])))
                mut_id = type + '|' + mut_ids
                mut_id_ls.append(mut_id)
    print('incorrect barcode, a total of %s reads were excluded' %error_read)
    print('deletion or insertion, a total of %s reads were excluded' %len_error)
    AA_dict = Counter(mut_id_ls)
    return AA_dict

def write_mut_table(mut_dic,outfilename):
  outfile = open(outfilename,'w')
  outfile.write("\t".join(['Mutation', 'count'])+"\n")
  for mut in mut_dic.keys():
    Mutation = mut
    count = mut_dic[mut]
    outfile.write("\t".join(map(str,[Mutation, count]))+"\n")
  outfile.close()
  print('Written %s' %outfilename)


def main():
    if len(sys.argv) != 4:
        sys.exit('[usage] python fastq2count.py <fastq file> < reference> <mutation table filename>')
    fastqfile = sys.argv[1]
    ref = sys.argv[2]
    outfilename = sys.argv[3]
    ref_dict = cal_ref_dict(ref)
    mutation_dic = cal_fastq_dic(fastqfile,ref_dict)
    write_mut_table(mutation_dic, outfilename)


if __name__ == "__main__":
  main()
