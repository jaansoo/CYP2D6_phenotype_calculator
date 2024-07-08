##BEFORE RUNNING:

#dependencies
import pandas as pd
import re
import math
import numpy as np
import argparse

parser=argparse.ArgumentParser(description="Process genotypes to determine CYP2D6 activity score and genotype inferred phenotype")
parser.add_argument("input_file",type=str,help="Path to the input TSV file grenerated by Cyrius tool by Illumina")
# parser.add_argument("separator",type=str,default='\t',help="the separator in your input file if using other type of input file")
parser.add_argument("output_file",type=str,help="Path to the putput CSV file")
parser.add_argument("star_allele_dict",type=str,help="path to the CSV file that contains the star allele to AS conversions")
args=parser.parse_args()

def mod_genotype(column):
    mod_list=[]
    for value in column:
        first_combination = value.split(';')[0]
        parts=first_combination.split('_')
        if len(parts)>2:
            first='+'.join(parts[:2])
            second='+'.join(parts[2:])
            result='/'.join([first,second])
            mod_list.append(result)
        elif len(parts)==2:
            result2='/'.join([parts[0],parts[1]])
            mod_list.append(result2)
        else:
            mod_list.append(value)
    return mod_list

def split_values_into_columns(column):
    Allele1=[]
    Allele2=[]
    for value in column:
        parts=value.split('/')
        if len(parts)<2:
            parts.append(None)
        Allele1.append(parts[0])
        Allele2.append(parts[1])
    return Allele1,Allele2

def get_activity_scores(star_allele,star_allele_dict):
    return star_allele_dict.get(star_allele,0)

def mapping_activity(allele1,allele2,star_allele_dict):
    #calculating the activity
    activity_list=[]
    for i in range(len(allele1)):
        activity_al1=0
        activity_al2=0
        #for allele1
        if 'x' in str(allele1[i]):
            a_value=str(allele1[i]).split("x")[0]
            times=str(allele1[i]).split("x")[1]
            times=int(times)
            a_parts=[a_value] * times
            a1='+'.join(a_parts)
        else:
            a1=allele1[i]
        
        if '+' in str(a1):
            b_value=str(a1).split("+")
            for val in b_value:
                acscore=get_activity_scores(val,star_allele_dict)
                activity_al1+=acscore
        else:
            activity_al1=get_activity_scores(str(a1),star_allele_dict)
        
        #for allele2
        if 'x' in str(allele2[i]):
            a_value=str(allele2[i]).split("x")[0]
            times=str(allele2[i]).split("x")[1]
            times=int(times)
            a_parts=[a_value] * times
            a2='+'.join(a_parts)
        else:
            a2=allele2[i]
        
        if '+' in str(a2):
            b_value=str(a2).split("+")
            for val in b_value:
                acscore=get_activity_scores(val,star_allele_dict)
                activity_al2+=acscore
        else:
            activity_al2=get_activity_scores(str(a2),star_allele_dict)
    
        activity_list.append([activity_al1,activity_al2])
    return pd.Series(activity_list)

def combining_alleles(allele1,allele2,star_allele_dict):
    activity_list=mapping_activity(allele1,allele2,star_allele_dict)
    activity=[]
    for item in activity_list:
        a,b=item
        if a is None or b is None or pd.isna(a) or pd.isna(b):
            activity.append('unknown')
        else:
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
        if item=='unknown':
            values.append('unknown')
        else:
            try:
                values.append(float(item))
            except ValueError:
                values.append('unknown')
    for value in values:
        if value == 'unknown':
            phenotype_list.append('unknown') #unknown type
        elif math.isclose(value,0.0,abs_tol=1e-9):
            phenotype_list.append('PM')# poor metabolisers
        elif 0 < value <1.25:
            phenotype_list.append('IM') # intermediate metabolisers
        elif 1.25 <= value <=2.25:
            phenotype_list.append('NM') #normal metabolisers
        elif value > 2.25:
            phenotype_list.append('UM') #ultrarapid metabolisers
        else:
            phenotype_list.append('unknown') 
    
    return pd.Series(phenotype_list)

#loading and reading the table file
df=pd.read_csv(args.input_file,sep='\t')

#modify genotype to include genotypes where no haplotypes has been asigned or which has more than one possible genotype
mod_genotypes=mod_genotype(df['Genotype'].tolist())
df['mod_genotype']=mod_genotypes

#split alleles if not already split
if 'Allele1' not in df.columns:
    Allele1,Allele2=split_values_into_columns(df['mod_genotype'])
    temp_df=pd.DataFrame({
        'Allele1':Allele1,
        'Allele2':Allele2
    })
    df[['Allele1','Allele2']]=temp_df

# Star_allele dictionary
dict_df=pd.read_csv(args.star_allele_dict)
star_allele_dict=dict_df.set_index('star_allele')['Activity_score'].to_dict()

#determining and writing the phenotype
df['Activity']=combining_alleles(df['Allele1'],df['Allele2'],star_allele_dict)
df['Cyrius_phenotype']=mapping_genotype_to_phenotype(df['Activity'])
df.to_csv(args.output_file,index=False)
print('Matching complete')