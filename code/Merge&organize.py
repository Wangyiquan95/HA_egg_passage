import os
import glob
import pandas as pd


files=glob.glob('result/*_count_singleAA.tsv')
print(files)
merge_df = pd.DataFrame()
for file in files:
    filename=(os.path.basename(file)).rsplit('_')[0]
    df = pd.read_table(file, sep='\t')
    df = df[df['count'] !=0]
    df['File_Name'] = filename
    merge_df=pd.concat([merge_df,df],axis=0)

info_df = pd.read_excel('data/sample information CPOS-211013-WL-9970a.xlsx',sheet_name=1)
DF=pd.merge(merge_df,info_df,on=['File_Name','Viral_RNA'],how='left')
DF['residues']=DF['Mutation'].str[1:-1]
DF['Group']=DF['Viral_RNA']+DF['Egg_passage']+DF['Biological_replicate']
DF['Freq']=((DF['count']) / DF.groupby(['Group','residues'])['count'].transform('sum'))
DF.to_excel('result/HA_EggPassaging_all_singleAAlevel.xlsx')