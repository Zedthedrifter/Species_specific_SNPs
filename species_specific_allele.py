#!/home/pholling/apps/env/snps/bin/python3
# 2025.08.12
# Zed ZChen@rbge.org.uk
# Input: vcf file, with GT:DP available (default)
# Output: csv file with genotype frequency for each species

from __future__ import division
import sys
import gzip
import copy
import argparse
import pandas as pd

def parse_vcf(vcf,output):
  #just in case the vcf file is gz
  if vcf.endswith(".gz"):
  	opener = gzip.open
  else:
  	opener = open
  #process each line
  snps,snps_info,n={},{},0
  for i in opener(vcf,'rt'):
    if i.startswith('#CHROM'):
      heading = i.strip().split()
      samples=heading[9:]
      samples=[s.split('/')[-1].replace('.bam','') for s in samples]
      print(f"processing {len(samples)} samples: {samples}")
    else:
      if not i.startswith('#') : #and n<=5000 : #/ this is a good cap for testing just a small number of SNPs
        n+=1
        snps[f"SNP_{n}"]={s:{} for s in samples}
        #COLLECT SNP INFO AND OUTPUT IN A CSV
        snps_info[f"SNP_{n}"]={i:j for i,j in zip(heading[:9],i.strip().split()[:9])}
        #COLLECT SAMPLE INFO
        i = i.strip().split()[9:] #the sample starts from teh 9th column
        genotypes=[gt.split(':')[0] for gt in i]
        prob=[gt.split(':')[1] for gt in i]
        for s,gt in zip (samples,genotypes):
          snps[f"SNP_{n}"][s]['GT']=gt
        for s,pl in zip (samples,prob):
          snps[f"SNP_{n}"][s]['PL']=pl.split(',')
  #OUTPUT SNP INFO: OKAY noboay wants to look at this aye?
  df=pd.DataFrame.from_dict(snps_info,orient='index') #turn to dataframe 
  df.to_csv(f'{output}/all_SNPs_summary.csv', index=True) 
  return(snps,snps_info,samples)

#get samples vs species
def parse_names(names):
  id_sp={l.split(',')[0].split('/')[-1].replace('.bam',''):l.split(',')[1].strip('\n') for l in open(names)}
  return(id_sp)

#calcualte GT per species
def species_specific_alleles(id_sp,snps,snps_info,output):
  sps=list(set(id_sp.values()))
  sp_nb=len(sps)
  snp_al,snp_GT,ssallele,hqallele,hqalGT={k:{} for k in snps.keys()},{k:{} for k in snps.keys()},{},{},{}
  for k,v in snps.items():
    sp_GT={}#a temporary dictionary to store genotype of each species, later we need to check if an allele is present in all sample if it can be species specific
    for sp in sps:
      gts=[v[sample]['GT'] for sample,spp in id_sp.items() if sp == spp]
      sp_GT[sp]=gts#store the genotype 
      alleles=[]
      for gt in gts:
        alleles+=[int(i) for i in gt.split('/') if i != '.'] #remove the uncertain ones
      #calculate allele frequency
      snp_al[k][sp]={al:round(alleles.count(al)/len(alleles)*100) for al in set(alleles)} #store the value
      
      #calculate GT frequency
      snp_GT[k][sp]={gt:round(100*gts.count(gt)/len(gts),2) for gt in set(gts)} #the dictionary with all GT frequencies for all SNPs
      
    tmp1={f"{species}_{al}": v.get(al,0) for species,v in snp_al[k].items() for al in [0,1]}
    tmp2={f"{species}_{gt}": v.get(gt,0) for species,v in snp_GT[k].items() for gt in ['0/0','0/1','1/1']}
    #print(snp_GT[k])
    for al,typ in {0:'REF',1:'ALT'}.items(): #assume two alleles
      freqs={sp:snp_al[k][sp].get(al,0) for sp in snp_al[k]}
      present=[1 for f in freqs.values() if f > 0 ]
      #first filter: if the allele is present in only one species:
      if sum(present)==1 : #the allele is present in only one species
        #get the species:
        target=[s for s,f in freqs.items() if f >0][0]
        ssallele[k]={}
        ssallele[k]['target sp']=target
        ssallele[k]['target allele']=snps_info[k][typ]
        ssallele[k]['target freq%']=freqs[target]
        ssallele[k]={**ssallele[k],**snps_info[k]}
        ssallele[k]={**ssallele[k],**tmp1}
        
        #second filter: if the allele is present in all samples of the species:
        present_in_all=[0 if str(al) in gt else 1 for gt in sp_GT[target]]
        if sum(present_in_all)==0: #if the allele is found in all samples, homo or hetero
          hqallele[k]=ssallele[k] #add allele freq info
          hqalGT[k]={**ssallele[k],**tmp2} #add GT freq info
          #print(hqallele[k])
          continue #no need to check the other allele
        else:
          continue
  #now select SNPs with species specific alleles
  df=pd.DataFrame.from_dict(ssallele,orient='index') #turn to dataframe 
  print(df)
  df.to_csv(f'{output}/specific_allele_freq.csv', index=True)  
  df=pd.DataFrame.from_dict(hqallele,orient='index') #turn to dataframe 
  print(df)
  df.to_csv(f'{output}/hq_specific_allele_freq.csv', index=True)  
  df=pd.DataFrame.from_dict(hqalGT,orient='index') #output allele as well as GT freq
  print(df)
  df.to_csv(f'{output}/hq_specific_allele_GT_freq.csv', index=True)  
      
        
      

#DONE WITH THE TRAINING DATASET
#ONTO THE QUERY SET
#EXTRACT GOOD SNPS
def extract_SNPs(snps,snps_info,snp_out):
  #reset the keys for extraction
  print('number of SNPs before renaming keys',len(snps),len(snps_info))
  snps={f"{v['#CHROM']}_{v['POS']}":snps[k] for k,v in snps_info.items()}
  snps_info={f"{v['#CHROM']}_{v['POS']}":v for k,v in snps_info.items()}
  print('number of SNPs after renaming keys',len(snps),len(snps_info))
  #extract
  print('number of species specific SNPs',len(snp_out)) 
  snps={k:snps.get(f"{v['#CHROM']}_{v['POS']}",'NA') for k,v in snp_out.items()}
  snps_info={k:snps_info.get(f"{v['#CHROM']}_{v['POS']}",'NA') for k,v in snp_out.items()}
  print('number of SNPs after extraction',len(snps),len(snps_info))
  return(snps,snps_info)

#SHOW SAMPLE GT
def specific_SNPs_GT(samples,snps,snps_info,snp_out,output):
  GTs={}
  for k,v in snp_out.items():
    if snps_info[k]=='NA': #if the snp is not called in those samples
      print(f'{k} not called in all samples: identical to reference')
    else:
      #print(snps_info[k],v)
      if snps_info[k]['ALT']==v['ALT']:
        GTs[k]={s:snps[k][s]['GT'] for s in samples}
      else:
        GTs[k]={s:snps[k][s]['GT'] if snps[k][s]['GT']=='0/0' else f"novo allele: {snps[k][s]['GT']}"  for s in samples}
#        GTs[k]={}
#        print(snps[k])
#        for s in snps[k]:
#          if snps[k][s]['GT']=='0/0':
#            GTs[k][s]='0/0'
#          else:
#            GTs[k]='novo'
  df=pd.DataFrame.from_dict(GTs,orient='index') #turn to dataframe 
  print(df)
  df.to_csv(f'{output}/specific_SNP_sample_GT.csv', index=True)  
  #convert GTs to list of alleles
  allele_table={0:'REF',1:'ALT'}
  alleles={snp:{} for snp in GTs}
  for snp,v in GTs.items():
    alleles[snp]={sample:[snps_info[snp][allele_table[int(i)]] for i in gt.split('/')] if 'novo' not in gt else f"novo allele:{[snps_info[snp][allele_table[int(i)]] for i in gt.replace('novo allele: ','').split('/')]}" for sample,gt in v.items()}
    alleles[snp]['target sp']=snp_out[snp]['target sp']
    alleles[snp]['target allele']=snp_out[snp]['target allele']
  df=pd.DataFrame.from_dict(alleles,orient='index') #turn to dataframe 
  print(df)
  df.to_csv(f'{output}/specific_SNP_sample_allele.csv', index=True)  
  return(GTs,alleles)

#TRANSLATE SAMPLE GT TO SPECIES
def GT_to_species(alleles,snp_out,output):
  species={}
  for k,v in alleles.items():
    species[k]={sample:snp_out[k]['target sp'] if snp_out[k]['target allele'] in al else 'NA' for sample,al in v.items()}
    species[k]['target sp']=snp_out[k]['target sp']
    species[k]['target allele']=snp_out[k]['target allele']
  df=pd.DataFrame.from_dict(species,orient='index') #turn to dataframe 
  print(df)
  df.to_csv(f'{output}/specific_allele_sample_species.csv', index=True)  
  return(species)

#PARENT FUNCTION 1
#using samples with known species, identify species specific SNPs
def specifi_SNPs(vcf,names,output):
  #print(vcf,names,output)
  snps,snps_info,samples=parse_vcf(vcf,output)
  id_sp=parse_names(names)
  snp_out=species_specific_alleles(id_sp,snps,snps_info,output) #output all species specific snps
  print('Complete calling species specific SNPs from training dataset')

#PARENT FUNCTION 2
#for samples without known species, find species/potential parents of hybrids based on species specific SNPs
def find_species(vcf,names,output):
  #read species specific SNPs from the previous step
  df=pd.read_csv(names,index_col=0)
  print(df)
  snp_out=df.to_dict(orient='index')
  #COLLECT SNPS
  snps,snps_info,samples=parse_vcf(vcf,output)
  #EXTRACT SPECIES SPECIFIC SNPS (IF EXIST)
  snps,snps_info=extract_SNPs(snps,snps_info,snp_out)
  GTs,alleles=specific_SNPs_GT(samples,snps,snps_info,snp_out,output)
  species=GT_to_species(alleles,snp_out,output)
  print('Complete calling species specific SNPs from query dataset and detecting parental strains')

#now making the execution function  
def main():
  count=0
  #parser = argparse.ArgumentParser(description='rename files of a certain type within the directory to use slurm array later',  usage = 'rename_files.py -i <indir> -o <outdir> --prefix <pre>')
  parser = argparse.ArgumentParser(description="Python script with multiple executable functions")
  subparsers = parser.add_subparsers(dest='command', help='Available commands')
  
  # Common arguments that will be reused
  common_args = argparse.ArgumentParser(add_help=False)
  common_args.add_argument('-v','--vcf', help='vcf file with GT and DP info',metavar='vcf')
  common_args.add_argument('-n','--names', help='sample ID to species names corresponding file, csv: sample ID/file path, species',metavar='names')
  common_args.add_argument('-o','--output',help='path to output directory',metavar='output')
  
  #find species specific SNPs
  func_parser=subparsers.add_parser('specifi_SNPs', parents=[common_args],help='find and output species specific SNPs', 
                                    usage = './calculate_snp.freq.py main_parent -v <vcf file> -n <ID_species.txt> -i <high bound> -l <low bound> -o <output_dir/prefix> ')
  func_parser.set_defaults(func=specifi_SNPs)
  
  #assign species to unclassified samples based on SNPs
  func_parser=subparsers.add_parser('find_species', parents=[common_args],help='find and output hits on identified species specific SNPs', 
                                    usage = './calculate_snp.freq.py find_species -v <vcf file> -n <species specific SNPs.csv> -i <high bound> -l <low bound> -o <output_dir/prefix> ')
  func_parser.set_defaults(func=find_species)
  
  #parse arguments
  args = parser.parse_args()

  if not args.command:
    parser.print_help()
    sys.exit(1)
  
  # Prepare common kwargs for the function call
  kwargs = {'vcf': args.vcf,
            'names': args.names,
            'output':args.output}
  
  # Add function-specific arguments
  #if args.command == 'rename_fqgz':
  #  kwargs['prefix'] = args.prefix
  #elif args.command == 'rename_fasta':
  #  kwargs['outfile'] = args.outfile
  #elif args.command == 'rename_contig':
  #  kwargs['infile'] = args.infile
  #  kwargs['fcsv'] = args.fcsv
  
  # Call the appropriate function
  args.func(**kwargs)
    

if __name__ == '__main__': 
  main()