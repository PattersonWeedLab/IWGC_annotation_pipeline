#!/bin/bash 
##==============================================
## files that get modified
genome="Eindica_glyres.genome.fa"
cupcake="E_ind_r.collapsed.gff"
##==============================================
## files that don't get modified
repeat_lib="${PWD}/Eleusine_indica_glyR-families.fa"
proteins="${PWD}/Ecoracana_560_v1.1.protein.fa"


##==============================================

## Split gff by chr
awk '{print >$1".gff"}' ${cupcake}

# chrs.list is list of tigs to be annotated.
for chr in `cat chrs.list` 
do 


## get name
module load SAMTools
if [ ! -d ${chr}  ] ; then mkdir ${chr} ; fi 
if [ ! -d "${PWD}/${chr}/tmp" ] ; then mkdir "${PWD}/${chr}/tmp" ; fi

## extract sequnce
fasta="${PWD}/${chr}/${chr}.fasta"

if [ -f "${genome}.fai" ]; then rm "${genome}.fai"; fi 
samtools faidx ${genome} ${chr} > ${fasta}
module unload SAMTools

## get transcripts
# Put Module load gffread here!
# OR create your own conda env that named as below and
source ~/.bashrc 
#conda activate iwgcfa001

module load gffread

rm "${genome}.fai"

transcripts="${PWD}/${chr}/${chr}_transcripts.fa"

gffread -g ${genome} -w ${transcripts} ${chr}.gff 


## Now add maker ctl scripts

## maker_bots.ctl

printf "\
#-----BLAST and Exonerate Statistics Thresholds
blast_type=ncbi+ #set to 'ncbi+', 'ncbi' or 'wublast'
use_rapsearch=0 #use rapsearch instead of blastx, 1 = yes, 0 = no

pcov_blastn=0.8 #Blastn Percent Coverage Threhold EST-Genome Alignments
pid_blastn=0.85 #Blastn Percent Identity Threshold EST-Genome Aligments
eval_blastn=1e-10 #Blastn eval cutoff
bit_blastn=40 #Blastn bit cutoff
depth_blastn=0 #Blastn depth cutoff (0 to disable cutoff)

pcov_blastx=0.5 #Blastx Percent Coverage Threhold Protein-Genome Alignments
pid_blastx=0.4 #Blastx Percent Identity Threshold Protein-Genome Aligments
eval_blastx=1e-06 #Blastx eval cutoff
bit_blastx=30 #Blastx bit cutoff
depth_blastx=0 #Blastx depth cutoff (0 to disable cutoff)

pcov_tblastx=0.8 #tBlastx Percent Coverage Threhold alt-EST-Genome Alignments
pid_tblastx=0.85 #tBlastx Percent Identity Threshold alt-EST-Genome Aligments
eval_tblastx=1e-10 #tBlastx eval cutoff
bit_tblastx=40 #tBlastx bit cutoff
depth_tblastx=0 #tBlastx depth cutoff (0 to disable cutoff)

pcov_rm_blastx=0.5 #Blastx Percent Coverage Threhold For Transposable Element Masking
pid_rm_blastx=0.4 #Blastx Percent Identity Threshold For Transposbale Element Masking
eval_rm_blastx=1e-06 #Blastx eval cutoff for transposable element masking
bit_rm_blastx=30 #Blastx bit cutoff for transposable element masking

ep_score_limit=20 #Exonerate protein percent of maximal score threshold
en_score_limit=20 #Exonerate nucleotide percent of maximal score threshold
"  >${chr}/maker_bopts.ctl


## maker_evm.ctl 
printf "\
#-----Transcript weights
evmtrans=10 #default weight for source unspecified est/alt_est alignments
evmtrans:blastn=0 #weight for blastn sourced alignments
evmtrans:est2genome=10 #weight for est2genome sourced alignments
evmtrans:tblastx=0 #weight for tblastx sourced alignments
evmtrans:cdna2genome=7 #weight for cdna2genome sourced alignments

#-----Protein weights
evmprot=10 #default weight for source unspecified protein alignments
evmprot:blastx=2 #weight for blastx sourced alignments
evmprot:protein2genome=10 #weight for protein2genome sourced alignments

#-----Abinitio Prediction weights
evmab=10 #default weight for source unspecified ab initio predictions
evmab:snap=10 #weight for snap sourced predictions
evmab:augustus=10 #weight for augustus sourced predictions
evmab:fgenesh=10 #weight for fgenesh sourced predictions
evmab:genemark=7 #weight for genemark sourced predictions " >${chr}/maker_evm.ctl

##maker_exe.ctl

printf "\
#-----BLAST and Exonerate Statistics Thresholds
blast_type=ncbi+ #set to 'ncbi+', 'ncbi' or 'wublast'
use_rapsearch=0 #use rapsearch instead of blastx, 1 = yes, 0 = no

pcov_blastn=0.8 #Blastn Percent Coverage Threhold EST-Genome Alignments
pid_blastn=0.85 #Blastn Percent Identity Threshold EST-Genome Aligments
eval_blastn=1e-10 #Blastn eval cutoff
bit_blastn=40 #Blastn bit cutoff
depth_blastn=0 #Blastn depth cutoff (0 to disable cutoff)

pcov_blastx=0.5 #Blastx Percent Coverage Threhold Protein-Genome Alignments
pid_blastx=0.4 #Blastx Percent Identity Threshold Protein-Genome Aligments
eval_blastx=1e-06 #Blastx eval cutoff
bit_blastx=30 #Blastx bit cutoff
depth_blastx=0 #Blastx depth cutoff (0 to disable cutoff)

pcov_tblastx=0.8 #tBlastx Percent Coverage Threhold alt-EST-Genome Alignments
pid_tblastx=0.85 #tBlastx Percent Identity Threshold alt-EST-Genome Aligments
eval_tblastx=1e-10 #tBlastx eval cutoff
bit_tblastx=40 #tBlastx bit cutoff
depth_tblastx=0 #tBlastx depth cutoff (0 to disable cutoff)

pcov_rm_blastx=0.5 #Blastx Percent Coverage Threhold For Transposable Element Masking
pid_rm_blastx=0.4 #Blastx Percent Identity Threshold For Transposbale Element Masking
eval_rm_blastx=1e-06 #Blastx eval cutoff for transposable element masking
bit_rm_blastx=30 #Blastx bit cutoff for transposable element masking

ep_score_limit=20 #Exonerate protein percent of maximal score threshold
en_score_limit=20 #Exonerate nucleotide percent of maximal score threshold
(iwgcfa001) [hallna12@watership makerp]$ cat maker_evm.ctl 
#-----Transcript weights
evmtrans=10 #default weight for source unspecified est/alt_est alignments
evmtrans:blastn=0 #weight for blastn sourced alignments
evmtrans:est2genome=10 #weight for est2genome sourced alignments
evmtrans:tblastx=0 #weight for tblastx sourced alignments
evmtrans:cdna2genome=7 #weight for cdna2genome sourced alignments

#-----Protein weights
evmprot=10 #default weight for source unspecified protein alignments
evmprot:blastx=2 #weight for blastx sourced alignments
evmprot:protein2genome=10 #weight for protein2genome sourced alignments

#-----Abinitio Prediction weights
evmab=10 #default weight for source unspecified ab initio predictions
evmab:snap=10 #weight for snap sourced predictions
evmab:augustus=10 #weight for augustus sourced predictions
evmab:fgenesh=10 #weight for fgenesh sourced predictions
evmab:genemark=7 #weight for genemark sourced predictions
(iwgcfa001) [hallna12@watership makerp]$ cat maker_exe.ctl
#-----Location of Executables Used by MAKER/EVALUATOR
makeblastdb=/opt/software/BLAST+/2.3.0/bin/makeblastdb #location of NCBI+ makeblastdb executable
blastn=/opt/software/BLAST+/2.3.0/bin/blastn #location of NCBI+ blastn executable
blastx=/opt/software/BLAST+/2.3.0/bin/blastx #location of NCBI+ blastx executable
tblastx=/opt/software/BLAST+/2.3.0/bin/tblastx #location of NCBI+ tblastx executable
formatdb=/opt/software/BLAST/2.2.26/bin/formatdb #location of NCBI formatdb executable
blastall=/opt/software/BLAST/2.2.26/bin/blastall #location of NCBI blastall executable
xdformat= #location of WUBLAST xdformat executable
blasta= #location of WUBLAST blasta executable
prerapsearch= #location of prerapsearch executable
rapsearch= #location of rapsearch executable
RepeatMasker=/opt/software/RepeatMasker/4.1.2/RepeatMasker #location of RepeatMasker executable
exonerate=/opt/software/exonerate/2.2.0/bin/exonerate

#-----Ab-initio Gene Prediction Algorithms
snap=/opt/software/SNAP/1.0/snap #location of snap executable
gmhmme3= #location of eukaryotic genemark executable
gmhmmp= #location of prokaryotic genemark executable
augustus= #location of augustus executable
fgenesh= #location of fgenesh executable
evm=/opt/software/EvidenceModeler/1.1.1/evidence_modeler.pl #location of EvidenceModeler executable
tRNAscan-SE=/opt/software/tRNAscan/2.0/tRNAscan-SE #location of trnascan executable
snoscan=/opt/software/snoscan/1.0/snoscan #location of snoscan executable

#-----Other Algorithms
probuild= #location of probuild executable (required for genemark)
" >${chr}/maker_exe.ctl 

## maker_opts.ctl

printf "\
#-----Genome (these are always required)
genome=${fasta} #genome sequence (fasta file or fasta embeded in GFF3 file)
organism_type=eukaryotic #eukaryotic or prokaryotic. Default is eukaryotic

#-----Re-annotation Using MAKER Derived GFF3
maker_gff= #MAKER derived GFF3 file
est_pass=0 #use ESTs in maker_gff: 1 = yes, 0 = no
altest_pass=0 #use alternate organism ESTs in maker_gff: 1 = yes, 0 = no
protein_pass=0 #use protein alignments in maker_gff: 1 = yes, 0 = no
rm_pass=0 #use repeats in maker_gff: 1 = yes, 0 = no
model_pass=0 #use gene models in maker_gff: 1 = yes, 0 = no
pred_pass=0 #use ab-initio predictions in maker_gff: 1 = yes, 0 = no
other_pass=0 #passthrough anyything else in maker_gff: 1 = yes, 0 = no

#-----EST Evidence (for best results provide a file for at least one)
est=${transcripts} #set of ESTs or assembled mRNA-seq in fasta format
altest= #EST/cDNA sequence file in fasta format from an alternate organism
est_gff= #aligned ESTs or mRNA-seq from an external GFF3 file
altest_gff= #aligned ESTs from a closly relate species in GFF3 format

#-----Protein Homology Evidence (for best results provide a file for at least one)
protein=${proteins}  #protein sequence file in fasta format (i.e. from mutiple organisms)
protein_gff=  #aligned protein homology evidence from an external GFF3 file

#-----Repeat Masking (leave values blank to skip repeat masking)
model_org= #select a model organism for RepBase masking in RepeatMasker
rmlib=${repeat_lib} #provide an organism specific repeat library in fasta format for RepeatMasker
repeat_protein=/opt/software/Maker/3.01.04/data/te_proteins.fasta #provide a fasta file of transposable element proteins for RepeatRunner
rm_gff= #pre-identified repeat elements from an external GFF3 file
prok_rm=0 #forces MAKER to repeatmask prokaryotes (no reason to change this), 1 = yes, 0 = no
softmask=1 #use soft-masking rather than hard-masking in BLAST (i.e. seg and dust filtering)

#-----Gene Prediction
snaphmm= #SNAP HMM file
gmhmm= #GeneMark HMM file
augustus_species= #Augustus gene prediction species model
fgenesh_par_file= #FGENESH parameter file
pred_gff= #ab-initio predictions from an external GFF3 file
model_gff= #annotated gene models from an external GFF3 file (annotation pass-through)
run_evm=0 #run EvidenceModeler, 1 = yes, 0 = no
est2genome=1 #infer gene predictions directly from ESTs, 1 = yes, 0 = no
protein2genome=1 #infer predictions from protein homology, 1 = yes, 0 = no
trna=1 #find tRNAs with tRNAscan, 1 = yes, 0 = no
snoscan_rrna= #rRNA file to have Snoscan find snoRNAs
snoscan_meth= #-O-methylation site fileto have Snoscan find snoRNAs
unmask=0 #also run ab-initio prediction programs on unmasked sequence, 1 = yes, 0 = no
allow_overlap= #allowed gene overlap fraction (value from 0 to 1, blank for default)

#-----Other Annotation Feature Types (features MAKER doesn't recognize)
other_gff= #extra features to pass-through to final MAKER generated GFF3 file

#-----External Application Behavior Options
alt_peptide=C #amino acid used to replace non-standard amino acids in BLAST databases
cpus=1 #max number of cpus to use in BLAST and RepeatMasker (not for MPI, leave 1 when using MPI)

#-----MAKER Behavior Options
max_dna_len=100000 #length for dividing up contigs into chunks (increases/decreases memory usage)
min_contig=1 #skip genome contigs below this length (under 10kb are often useless)

pred_flank=200 #flank for extending evidence clusters sent to gene predictors
pred_stats=0 #report AED and QI statistics for all predictions as well as models
AED_threshold=1 #Maximum Annotation Edit Distance allowed (bound by 0 and 1)
min_protein=0 #require at least this many amino acids in predicted proteins
alt_splice=0 #Take extra steps to try and find alternative splicing, 1 = yes, 0 = no
always_complete=0 #extra steps to force start and stop codons, 1 = yes, 0 = no
map_forward=0 #map names and attributes forward from old GFF3 genes, 1 = yes, 0 = no
keep_preds=0 #Concordance threshold to add unsupported gene prediction (bound by 0 and 1)

split_hit=10000 #length for the splitting of hits (expected max intron size for evidence alignments)
min_intron=20 #minimum intron length (used for alignment polishing)
single_exon=0 #consider single exon EST evidence when generating annotations, 1 = yes, 0 = no
single_length=250 #min length required for single exon ESTs if 'single_exon is enabled'
correct_est_fusion=0 #limits use of ESTs in annotation to avoid fusion genes

tries=2 #number of times to try a contig if there is a failure for some reason
clean_try=0 #remove all data from previous run before retrying, 1 = yes, 0 = no
clean_up=0 #removes theVoid directory with individual analysis files, 1 = yes, 0 = no
TMP="${PWD}/${chr}/tmp" #specify a directory other than the system default temporary directory for temporary files
" >${chr}/maker_opts.ctl

done 



