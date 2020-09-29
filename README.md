# SCcaller_NB
Single Cell Caller (SCcaller) - NB fork - Identify single nucleotide variations (SNVs) and short insertions and deletions (INDELs) from single cell sequencing data

Version 2.0.0_NB

Updated date: 2020.09.24

Cite us:

Dong X et al. Accurate identification of single-nucleotide variants in whole-genome-amplified single cells. Nat Methods. 2017 May;14(5):491-493. doi: 10.1038/nmeth.4227.

#####
## Author and License

Authors: Xiao Dong, Yujue Wang

Email: biosinodx@gmail.com (X.D.), xiao.dong@einstein.yu.edu (X.D.), spsc83@gmail.com (Y.W.)

Modifications: Nico Borgsmüller (nico.borgsmueller@bsse.ethz.ch)

Licensed under the GNU Affero General Public License version 3 or later

#####
## DEPENDENCIES

python 2.7 or 3.X (not working with vcf catalog and bed output)

numpy
pysam(0.15.1) or samtools v.1.9+ (Other versions not tested)


#####
## FORK MODIFICATIONS

* Adaptation for python 3.X
* Implementation redone, simplified, and variables renamed for easier maintanance 
* New arguments added: --minvarfrac, --debugging, --seed
* If 'multiple-genotype' filter is set, all fulfill variant reads > minvar reads
* FORMAT/GT is 0/1 now, if likelihood for it is highest and #var reads > #total reads * minvarfrac (previously: highest likelihood and 7 * #var reads > #ref reads)
* Likelihood for wildtype adjusted (1/8 * theta instead of 1/8)
* FORMAT/BN added to provide more detail why SO is True|False|NA



#####
## USAGE

###
### I. Basic usage: calling SNVs and INDELs from a cell

#### I.a When you have heterozygous SNPs pre-called from bulk DNA of the same subject,

python sccaller_v2.0.0.py \

  --bam cell.bam \ # bam file of a single cell
  
  --fasta ref.fa \ # reference genome in fasta format
  
  --output cell.vcf \ # output vcf file
  
  --snp_type hsnp \ # using heterozygous SNPs pre-called from bulk DNA
  
  --snp_in hsnp.vcf (or bed) \ # vcf or bed file of heterozygous SNPs pre-called from bulk DNA
  
  --cpu_num 8 \ # using 8 cpu threads
  
  --engine samtools # using samtools engine

#### I.b When you do not have heterozygous SNPs pre-called from bulk DNA of the same subject, obtain SNPs from dbSNP or other databases,

python sccaller_v2.0.0.py \

  --bam cell.bam \ # bam file of a single cell
  
  --fasta ref.fa \ # reference genome in fasta format
  
  --output cell.vcf \ # output vcf file
  
  --snp_type dbsnp \ # using SNPs from dbSNP database (or other database)
  
  --snp_in dbsnp.vcf (or bed) \ # vcf or bed file containing all SNPs in dbSNP (or other) database
  
  --cpu_num 8 \ # using 8 cpu threads
       
  --engine samtools # using samtools engine

### II. Calling somatic SNVs and INDELs not present in bulk DNA

#### II.a Step 1. Calling SNVs and INDELs from a cell together with bulk DNA in input,

python sccaller_v2.0.0.py \

  --bam cell.bam \ # bam file of a single cell
  
  --bulk bulk.bam \ # bam file of bulk DNA
  
  --fasta ref.fa \ # reference genome in fasta format
  
  --output cell.vcf \ # output vcf file
  
  --snp_type hsnp \ # using heterozygous SNPs pre-called from bulk DNA
  
  --snp_in hsnp.vcf (or bed) \ # vcf or bed file of heterozygous SNPs pre-called from bulk DNA
  
  --cpu_num 8 \ # using 8 cpu threads
     
  --engine samtools # using samtools engine

#### II.b Step 2. Filtering out SNVs and INDELs observed in bulk or sequencing depth <= 20x in the single cell

grep '0/1' cell.vcf | grep 'True' | awk '$7=="." && length($5)==1' | awk -F "[:,]" '$8+$9>=20' > cell.somatic.snv.vcf

grep '0/1' cell.vcf | grep 'True' | awk '$7=="." && length($5)>1' | awk -F "[:,]" '$8+$9>=20' > cell.somatic.indel.vcf

### III. Notes on X/Y chromosome in males and ploidy 

Please note, sccaller was designed assuming diploidy genome and two copies of chromosomes. It cannot be used for calling mutations from X/Y chromosome of a male subject.
