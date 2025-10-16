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
      if not i.startswith('#'): #and n<=100; this is a good cap for testing just a small number of SNPs
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
def GT_per_species2(id_sp,snps,snps_info,high,low,output):
  log=open(f"{output}/log.txt",'w')
  log.write(f"species specific SNP with a GT\n>= {high}% frequency in the present species,and \n<{low}% in absent species\n")
  sps=list(set(id_sp.values()))
  log.write(f"species: {sps}\n")
  sp_nb=len(sps)
  GTs=set([v2['GT'] for v in snps.values() for v2 in v.values()]) #extract all possible genotypes
  GTs={gt: f"GT_{i}" for i,gt in enumerate(GTs)}
  print(f"genotypes: {GTs}")
  log.write(f"genotypes: {GTs}\n")
  snp_GT={k:{} for k in snps.keys()}
  print(f"species included in this table: {sps}")
  for k,v in snps.items():
    for sp in sps:
      gts=[v[sample]['GT'] for sample,spp in id_sp.items() if sp == spp]
      snp_GT[k][sp]={gt:round(100*gts.count(gt)/len(gts),2) for gt in set(gts)} #the dictionary with all GT frequencies for all SNPs
  GT_dict={gt:{} for gt in GTs}
  #look at each case of the GT: 0/0,0/1,1/1,etc.
  snp_hits={}
  for gt,name in GTs.items():
    gt_freq={k:{species:freq.get(gt,0) for species,freq in v.items()} for k,v in snp_GT.items()} 
    log.write(f"process {len(gt_freq)} SNPs\n")
    GT_dict[gt]=gt_freq #the frequency of each genotype of each species
    df=pd.DataFrame.from_dict(gt_freq,orient='index') #turn to dataframe 
    df.to_csv(f'{output}/{name}_freq.csv', index=True) 
    #print(df)
    #SELECT FOR SPECIES SPECIFIC GT OF THE SNP
    selected={}
    for k,v in gt_freq.items():
      freqs=[f for sp,f in v.items()]
      over=[1 for i in freqs if i >= high]
      belo=[1 for i in freqs if i <=  low ]
      if sum(over)==1 and sum(over)+sum(belo)==sp_nb: #only one species has 'over' and 'over'+'below' adds up to the species number
        snp_hits[k]=gt
        selected[k]=v
        #add SNP info 
        selected[k]['genotype']=gt
        selected[k]={**selected[k],**snps_info[k]} #python3.5 or greater
    log.write(f"found {len(selected)} SNPs with species specific GT {gt}\n")
    df=pd.DataFrame.from_dict(selected,orient='index') #turn to dataframe 
    print(df)
    df.to_csv(f'{output}/{name}_freq_{low}_{high}.csv', index=True)  
  log.close() 
  #get all species specific SNPs
  snp_out={snp:{f"{species}_{gt}": snp_GT[snp][species].get(gt,0) for species,v in snp_GT[snp].items() for gt in GTs} for snp in snp_hits}
  species=set(id_sp.values())
  for k in snp_out:
    target=[sp for sp in species if snp_out[k][f"{sp}_{snp_hits[k]}"]>=high]
    snp_out[k]['species specific GT']=snp_hits[k]
    snp_out[k]['target species']=target[0]
    snp_out[k]={**snp_out[k],**snps_info[k]}
  df=pd.DataFrame.from_dict(snp_out,orient='index') #turn to dataframe 
  print(df)
  df.to_csv(f'{output}/specific_SNP_freq_{low}_{high}.csv', index=True)  
  return(snp_out)

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
  print(snps)
  return(snps,snps_info)

#SHOW SAMPLE GT
def specific_SNPs_GT(samples,snps,snps_info,snp_out,output):
  GTs={}
  for k,v in snp_out.items():
    if v['species specific GT']!='./.':
      if snps_info[k]=='NA': #if the snp is not called in those samples
        #GTs[k]={s:'0/0' for s in samples}
        print(f'{k} not called in all samples: identical to reference')
      else:
        #print(snps_info[k],v)
        if snps_info[k]['ALT']==v['ALT']:
          GTs[k]={s:snps[k][s]['GT'] for s in samples}
        else:
          for s in snps[k]:
            if snps[k][s]['GT']=='0/0':
              GTs[k][s]='0/0'
            else:
              GTs[k]='novo'
  df=pd.DataFrame.from_dict(GTs,orient='index') #turn to dataframe 
  print(df)
  df.to_csv(f'{output}/specific_SNP_sample_GT.csv', index=True)  
  return(GTs)

#TRANSLATE SAMPLE GT TO SPECIES
def GT_to_species(GTs,snp_out,output):
  species={}
  for k,v in GTs.items():
    species[k]={sample:snp_out[k]['target species'] if gt==snp_out[k]['species specific GT'] else 'NA' for sample,gt in v.items()}
    species[k]['target species']=snp_out[k]['target species']
    species[k]['species specific GT']=snp_out[k]['species specific GT']
  df=pd.DataFrame.from_dict(species,orient='index') #turn to dataframe 
  print(df)
  df.to_csv(f'{output}/specific_SNP_sample_species.csv', index=True)  
  return(species)

#PARENT FUNCTION 1
#using samples with known species, identify species specific SNPs
def specifi_SNPs(vcf,names,high,low,output):
  #print(vcf,names,output)
  snps,snps_info,samples=parse_vcf(vcf,output)
  id_sp=parse_names(names)
  snp_out=GT_per_species2(id_sp,snps,snps_info,int(high),int(low),output) #output all species specific snps

#PARENT FUNCTION 2
#for samples without known species, find species/potential parents of hybrids based on species specific SNPs
def find_species(vcf,names,high,low,output):
  #read species specific SNPs from the previous step
  df=pd.read_csv(names,index_col=0)
  print(df)
  snp_out=df.to_dict(orient='index')
  #COLLECT SNPS
  snps,snps_info,samples=parse_vcf(vcf,output)
  #EXTRACT SPECIES SPECIFIC SNPS (IF EXIST)
  snps,snps_info=extract_SNPs(snps,snps_info,snp_out)
  GTs=specific_SNPs_GT(samples,snps,snps_info,snp_out,output)
  species=GT_to_species(GTs,snp_out,output)
  print('done')

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
  common_args.add_argument('-i','--high', help='GT freq of species specific GT to be at least xx% in the present species',metavar='high')
  common_args.add_argument('-l','--low', help='GT freq of species specific GT to be at most xx% in the absent species',metavar='low')
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
            'high':args.high,
            'low':args.low,
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