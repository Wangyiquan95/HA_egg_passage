import pandas as pd
from Bio import SeqIO,bgzf
import gzip
import glob
from collections import defaultdict


def Call_barcode_type(metadata_file):
    fq_df = pd.read_excel(metadata_file)
    rename_dict=defaultdict(dict)
    for index, info_dict in fq_df.iterrows():
        fq1name=info_dict['original_filename1']
        fq2name = info_dict['original_filename2']
        barcode=info_dict['3NT_BARCODE']
        # write fqfilename to dict[orignal filename][barcode]:new filenames
        rename_dict[fq1name][barcode] = info_dict['renamed_filename1']
        rename_dict[fq2name][barcode] = info_dict['renamed_filename2']
    return rename_dict

def demultiplex_fastq(file_path,rename_dict):
    demulti_dict = defaultdict(list)
    print ("reading %s" % file_path)
    handle_in = gzip.open(file_path, "rt")
    Rrecords = SeqIO.parse(handle_in,"fastq")
    fastq=file_path.split('/')[1]
    for record in Rrecords:
        seq = str(record.seq)
        #filter out ambiguous read
        # if 'N' in seq:continue
        barcode = seq[:3]
        if barcode in rename_dict[fastq].keys():
            new_name = rename_dict[fastq][barcode]
            demulti_dict[new_name].append(record)
    return demulti_dict

def write_fastq(demulti_dict):
    for fastq in demulti_dict.keys():
        out_file='cleaned_fastq/'+fastq
        with bgzf.BgzfWriter(out_file, "wb") as output_handle:
            SeqIO.write(sequences=demulti_dict[fastq], handle=output_handle, format="fastq")
    print('Written %s' %out_file)

def main():
    metafile='data/sra metadata_v4.xlsx'
    rename_dict = Call_barcode_type(metafile)
    files = glob.glob('fastq/*fastq.gz')
    for file in files:
        demulti_dict=demultiplex_fastq(file,rename_dict)
        write_fastq(demulti_dict)

if __name__ == "__main__":
  main()
