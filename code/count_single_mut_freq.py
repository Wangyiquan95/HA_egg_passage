import os
import glob
import pandas as pd

def cal_mut_ls(df):
    type_ls = list(set(df.Viral_RNA.tolist()))
    mut_set = []
    muts = df.Mutation.tolist()
    for mut in muts:
        if '-' in mut:
            a_mut = mut.split('-')
            mut_set.extend(a_mut)
        else:
            mut_set.append(mut)
    mut_ls = list(set(mut_set))
    return mut_ls,type_ls
files=glob.glob('result/*_count.tsv')

for file in files:
    filename=os.path.basename(file).split('.')[0]
    df = pd.read_table(file, sep='\t')
    df = df[df['count'] >= 10]
    df[['Viral_RNA','Mutation']] = df['Mutation'].str.split('|',expand=True)
    mut_ls, type_ls = cal_mut_ls(df)
    output='result/'+filename+'_singleAA.tsv'
    print('working on ',filename)
    with open (output,'w') as f:
        header=['Viral_RNA','Mutation','count']
        f.write('\t'.join(header)+'\n')
        for type in type_ls:
            for mut in mut_ls:
                mut_df = df[(df.Viral_RNA == type)&(df.Mutation.str.contains(mut))]
                mut_count = str(mut_df['count'].sum())
                f.write('\t'.join([type,mut,mut_count])+'\n')
