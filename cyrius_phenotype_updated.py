##BEFORE RUNNING:
## change the data file to use OR comment in/out the files as needed
## change the output file to match the necessary file

#loading dependencies
import pandas as pd
import re
import math
import sys
import numpy as np

print(f"List of arguments: {sys.argv}")

#loading and reading the table file
df=pd.read_csv(sys.argv[1],sep='\t')  # for bb data
# print(df,df.columns)
# df = pd.read_table('cyrius genotyping out/cyrius_output.tsv',sep='\t') # for 50genomes data

# def remove_(column):
#     mod_list=[]
#     for value in column:
#         parts=value.split('_')
#         if len(parts)>2:
#             first='+'.join(parts[:2])
#             second='+'.join(parts[2:])
#             result='/'.join([first,second])
#             mod_list.append(result)
#         else:
#             mod_list.append(value)

# mod_genotype=remove_(df['Genotype'].tolist())
# df['mod_genotype']=mod_genotype


def split_values_into_columns(column):
    Allele1=[]
    Allele2=[]
    for value in column:
        first_combination = value.split(';')[0]
        parts=first_combination.split('/')
        if len(parts)<2:
            parts.append(None)
        Allele1.append(parts[0])
        Allele2.append(parts[1])
    return Allele1,Allele2


#split alleles if not already split
if 'Allele1' not in df.columns:
    Allele1,Allele2=split_values_into_columns(df['Genotype'])
    temp_df=pd.DataFrame({
        'Allele1':Allele1,
        'Allele2':Allele2
    })
    df[['Allele1','Allele2']]=temp_df

#defining the functions
#TODO  add more star alleles and their values
functional_genotypes = ['*1','*2','*2+*4','*13+*2','*27','*35','*34','*33','*39','*45'] #1.0
reduced_functionality_alleles = ['*17','*29','*41x2','*59','*10+*9','*10x2'] #0.5
severely_reduced_alleles=['*36+*10','*10','*32','*41','*9','*10+*4','*119'] #0.25                              
non_functional_alleles = ['*3','*4','*5','*6','*6x2','*7','*11','*13','*15','*40','*68','*68+*4','*68+*68','*68+*68+*4',
                            '*4.013+*4','*4.013+*4.013','*4.013','*4x2','*68+*68+*68+*68+*4'] #0
unknown_undetermined_alleles=['*22','*24','*28','*43','*43x2','*83','*106','*108','*116','*117','*122','*127','*131']   #need to add a calculation
increased_function_alleles= ['*1x2','*2x2','*1+*1','*35x2','*33x2'] #2.0
x3_increase=['*2x3']#3.0
x4_increase=['*2x4']#4.0
def mapping_activity(allele1,allele2):
    #calculating the activity
    activity_list=[]
    for i in range(len(allele1)):
        activity_al1=0
        activity_al2=0
        if str(allele1[i]) in functional_genotypes:
            activity_al1='1.00'
        elif str(allele1[i]) in reduced_functionality_alleles:
            activity_al1='0.50'
        elif str(allele1[i]) in severely_reduced_alleles:
            activity_al1='0.25'
        elif str(allele1[i]) in non_functional_alleles:
            activity_al1='0.00'
        elif str(allele1[i]) in unknown_undetermined_alleles:
            activity_al1='*'
        elif str(allele1[i]) in increased_function_alleles:
            activity_al1='2.00'
        elif str(allele1[i]) in x3_increase:
            activity_al1='3.00'
        elif str(allele1[i]) in x4_increase:
            activity_al1='4.00'
        else:
            activity_al1='-1'
        if str(allele2[i]) in functional_genotypes:
            activity_al2='1.00'
        elif str(allele2[i]) in reduced_functionality_alleles:
            activity_al2='0.50'
        elif str(allele2[i]) in severely_reduced_alleles:
            activity_al2='0.25'
        elif str(allele2[i]) in non_functional_alleles:
            activity_al2='0.00'
        elif str(allele2[i]) in unknown_undetermined_alleles:
            activity_al2='*'
        elif str(allele2[i]) in increased_function_alleles:
            activity_al2='2.00'
        elif str(allele2[i]) in x3_increase:
            activity_al2='3.00'
        elif str(allele2[i]) in x4_increase:
            activity_al2='4.00'
        else:
            activity_al2='-1'
        activity_list.append([activity_al1,activity_al2])
    return pd.Series(activity_list)

def combining_alleles(allele1,allele2):
    activity_list=mapping_activity(allele1,allele2)
    activity=[]
    for item in activity_list:
        a,b=item
        try:
            activity.append(str(float(a) +float(b)))
        except ValueError:
            activity.append(str(a+b))
    return activity


def mapping_genotype_to_phenotype(activity_list):
    #assigning the phenotype
    phenotype_list=[]
    values=[]
    for item in activity_list:
        matches = re.findall(r'(-?\d+\.?\d*)\*?', item)
        if matches:
            values.append(float(matches[0]))
        else:
            values.append(np.nan)
    for value in values:
        if value >=2.25:
            phenotype_list.append('UM') #ultrafast metabolisers
        elif 1 <= value  < 2.25:
            print('adding normal value')
            phenotype_list.append('NM') #normal metabolisers
        elif 0 < value <1:
            phenotype_list.append('IM') # intermediate metabolisers
        elif math.isclose(value,0.0,abs_tol=1e-9):
            phenotype_list.append('PM')# poor metabolisers
        elif value =='aa':
            phenotype_list.append('missing') 
        else:
            phenotype_list.append('missing') 
    
    return pd.Series(phenotype_list)

#determining and writing the phenotype
df['Activity']=combining_alleles(df['Allele1'],df['Allele2'])
df['Cyrius_phenotype']=mapping_genotype_to_phenotype(df['Activity'])
df.to_csv(sys.argv[2],index=False)
print('Matching complete')