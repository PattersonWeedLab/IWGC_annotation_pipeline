# Current annotation pipeline of the International Weed Genomics Consortium ([IWGC](https://www.weedgenomics.org/))  
**Authors: [Dr. Nathan D. Hall](https://github.com/NDHall) (Developer), [Nicholas A. Johnson](https://github.com/Scrumpis) (Documentation), and [Dr. Eric L. Patterson](https://github.com/PattersonWeedLab) (P.I.)**  
**GitHub Repo: [IWGC_annotation_pipeline](https://github.com/PattersonWeedLab/IWGC_annotation_pipeline)**  



# Introduction

## About
This repostiory was established to document the genome annotation methods used by the [International Weed Genomics Consortium](https://www.weedgenomics.org/) and in the article "Subtelomeric 5-enolpyruvylshikimate-3-phosphate synthase copy number variation confers glyphosate resistance in Eleusine indica" and to also provide a publicly available pipeline or framework for computational genome annotation. The current state of the pipeline works well for the [Patterson Lab](https://github.com/PattersonWeedLab) but will require a substantial time investment for tweaking of small things like directory paths of files or programs and subsequent troubleshooting during the initial setup. It should also be noted that for functional annotation the [MetaCyc](https://metacyc.org/) database require a usage license (intended for academic use) that will be required for you to obtain or you will have to omit these parts to use the pipeline. We hope at the very least this documentation will provide researchers with a reference for designing their own genome annotation pipeline. We may attempt to release a more distributable version in the future, but there are no set plans currently.

  
## Summary of Structural Annotation

### 1. Handling Repeats
* Repeat regions are annotated using [RepeatModeler](https://github.com/Dfam-consortium/RepeatModeler) (for program versions see [Dependencies](#Dependencies)).
  
* Annotated repeat regions are masked with [RepeatMasker](https://github.com/rmhubley/RepeatMasker) to reduce
  the computational burden of further analysis.

* [bedtools](https://github.com/arq5x/bedtools2) is used to soft mask the genome using the output
  of RepeatMasker.
  
### 2. Map Isoseq Reads to Masked Genome
#### minimap2 alignment
* IsoSeq reads are mapped to the repeat-masked genome using [minimap2](https://github.com/lh3/minimap2).

#### SAMtools SAM to BAM
* [SAMtools](https://github.com/samtools/samtools) is used to convert the minimap2 alignment SAM file into a BAM file.

### 3. Collapse Isoforms
#### Cupcake:
* Isoforms are collapsed with [Cupcake](https://github.com/Magdoll/cDNA_Cupcake).  

### 4. Maker
* Genome, collapsed cDNA from Cupcake, repeat libraries from RepeatModeler, and a protein FASTA from a close relative species are fed into [MAKER](https://github.com/Yandell-Lab/maker).

### 5. Merge and Cleanup
* Genes that produce proteins under 27 amino acids long are removed from further annotation with only the longest proteins from each gene and unique untranslated regions (UTRs) used for functional annotation.


## Summary of Functional Annotation

### 1. Isoform Selection
* Isoforms are selected using [AGAT](https://github.com/NBISweden/AGAT) which 
selects the longest isoforms with longest intact open read frame.
  
* Isoforms are then extracted from [GFF3](https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md)
  format using using [gffread](https://github.com/gpertea/gffread)
  
* These predicted proteins are used for all downstream analysis.

### 2. Location Prediction
* [MultiLoc2 Workstation Edition](https://github.com/NDHall/MultiLoc2-1/tree/workstation-edition) (ML2) is used predict protein location within the cell. Briefly, ML2 uses a trained machine learning classifer to predict protein location. For more on this program see the excellent work by [Blum et al. (2009)](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-10-274).  

### 3. InterPro Analysis
* InterPro is run locally using iprscan5 version 5.47-82.0.
  
* It uses several programs to predict InterPro Accessions which can be broken down into the following categories.
    
    1. Sites
        * Conserved Site
        * Active Site
        * Binding Site
        * PTM
    2. Repeats
    3. Domains
    4. Families
    5. Homologous SuperFamilies
       
* For more information see release notes [here](https://www.ebi.ac.uk/interpro/release_notes/84.0/).
  
* IPRSCAN Results provide [MetaCyc](https://metacyc.org/) accessions but not descriptions. To obtain descriptions a custom python script was used to extract Pathway IDs and link them to their descriptions. MetaCyc is a proprietary database that is made conditionally available to academic researchers free of charge.  
  
### 4. Homology Search
*This search is for the best hit in each database using MMSeqs2. Which is an extremley fast method of searching that can be up to 36 times faster than BLAST [(Steinegger and Söding 2017)](https://www.nature.com/articles/nbt.3988).  

#### UniRef50
* UniRef50 is a clustered database that is regularly updated. See work by [Suzek et al. (2015)](https://academic.oup.com/bioinformatics/article/31/6/926/214968).
* We use [UniRef_50 version 2021_03](https://ftp.uniprot.org/pub/databases/uniprot/previous_releases/release-2021_03/relnotes.txt) and downloaded representative plants from the same release to add descriptions.  
Cluster information and UniRef IDs are reported

#### NCBI
All available proteins from fully sequenced plants were downloaded
May 2021 and put into a MMeqs2 database using [NCBI Datasets](https://www.ncbi.nlm.nih.gov/datasets/) and MMSeqs2 tool

### 5. Herbicide Interacting Genes
Herbicide interacting genes were identified from the literature and mapped to KEGG Ortholog. UniRef50 IDs
were converted to KEGG Orthologs and cross matched with herbicide interacting genes. Matches were flagged and 
Annotated with HRAC and WSSA classifiers.

### Output

The pipeline runs all analysis seperately, but a series of python scripts is used to collect, filter and 
and format the data. Results are reported in JSON format which is ideal for NOSQL databases and queries.
Given the number of Null values or missing values this format is ideal. The results are gappy because not every analysis 
returns a result for each isoform. This is expected. Finally, GFF3 is produced that contains Ontology Terms, and
Notes.

[fullsize poster](https://github.com/NDHall/NCWSS_2021/blob/main/media/NCWSS_poster_2021.pdf)

#### Small Poster for Reference

![thumbnail](https://github.com/NDHall/NCWSS_2021/blob/main/media/NCWSS_poster.png)

----------------


# Installation

## Dependencies
The below program version numbers are the exact versions used by the current IWGC_annotation_pipeline. Different versions may also be compatible but no other versions have been verified to be compatible with this pipeline by the PattersonWeedLab.

### Structural Annotation Dependencies:
* [RepeatModeler](https://github.com/Dfam-consortium/RepeatModeler) (version 2.0.2)
* [h5py](https://github.com/h5py/h5py) (version )
* [RepeatMasker](https://github.com/rmhubley/RepeatMasker) (version 4.1.2)
* [bedtools](https://github.com/arq5x/bedtools2) (version 2.30.0)
* [minimap2](https://github.com/lh3/minimap2) (version 2.24)
* [SAMtools](https://github.com/samtools/samtools) (version 1.11)
* [Cupcake](https://github.com/Magdoll/cDNA_Cupcake) (version 28.0)
* [gffread](https://github.com/gpertea/gffread) (version 0.12.7)
* [MAKER](https://github.com/Yandell-Lab/maker) (version 3.01.04)
* [gff3](https://pypi.org/project/gff3/) (version 1.0.1)
* [BCBio](https://pypi.org/project/bcbio-gff/) (version 0.7.0)
* [gfftools](https://github.com/ihh/gfftools) (version )

### Functional Annotation Dependencies:
* [AGAT](https://github.com/NBISweden/AGAT) (version 0.8.0)
* [gffread](https://github.com/gpertea/gffread) (version 0.12.7)
* [InterProScan5](https://github.com/ebi-pf-team/interproscan) (version 5.47-82.0)
* [MMseqs2](https://github.com/soedinglab/MMseqs2) (version 4.1)
* [MultiLoc2 Workstation Edition](https://github.com/NDHall/MultiLoc2-1/tree/workstation-edition) (version 2-1)

## Functional Annotation Scripts:
`git clone `

----------------


# Usage
Genomes are first structurally annotated in terms of repeat regions and gene models. Structural annotation is then used to functionally annotate gene models.

## Structural Annotation Usage

### 1. Handle Repeats
#### RepeatModeler:
First build a database for RepeatModeler.  
```bash
BuildDatabase -name Species_Name genome.fasta
``` 

Next use repeat database from last step to model repeats. -pa 1 uses 4 threads, -pa 25 uses 100 threads. Adjust accordingly. 
```bash
RepeatModeler -database Species_Name -pa 25 -LTRStruct
```

#### RepeatMasker:
Output from RepeatModeler used as input for RepeatMasker.  
```bash
RepeatMasker -gff -a -pa 20 -u -lib RepeatModeler_output-families.fa genome.fasta
```  
May need to load module h5py if the above fails.

#### bedtools:
Output from RepeatMasker used to soft mask repeat regions in genome with bedtools.  
```bash
bedtools maskfasta -fi genome.fasta -bed RepeatMasker_output.gff -soft -fo genome.softmasked.fasta
```

### 2, Map Isoseq Reads to Masked Genome
#### minimap2 alignment  
```bash
minimap2 -a -x splice -H -t 100 -O6,24 -B4 genome.softmasked.fasta isoseq.fastq -o output.sam
```

#### SAMtools
Prepare algnment for CupCake with samtools.
```bash
samtools view -b -T original_genome.fasta minimap2_alignment.sam > minimap2_alignment.bam
```

### 3. Collapse Isoforms with Cupcake
#### CupCake:
Collapse isoforms with CupCake.
```bash
python path/to/collapse_isoforms_by_sam.py --input ISOseq.fq --fq -b minimap2.sorted.bam -o Species_Name
```

### 4. MAKER (parallel run)
This is our workaround to parallelize MAKER.
#### Create a list of all FASTA headers to loop through
```bash
grep '>' genome.fasta | sed 's/>//g' > chrs.list
```

#### maker_split_and_run.sh
For intital setup, ensure all paths to MAKER required executables are accurate on your system in the below section.  
`#-----Location of Executables Used by MAKER/EVALUATOR`
  
Edit maker_split_and_run.sh for your input files as shown below.
```bash
#!/bin/bash 
##==============================================
## files that get modified
genome="Species_Name.softmasked.fasta"
cupcake="Species_Name.collapsed.gff"
##==============================================
## files that don't get modified
repeat_lib="Species_Name-families.fa"
proteins="Related_Species_Name.proteins.fa"
```

Run maker_split_and_run.sh to separate the genome into different directories by chromosome/scaffold/contig/number of headers in fasta and setup a MAKER.ctl for each.  
```bash
bash maker_split_and_run.sh
```

#### maker_run.sh
Edit the path `/path/to/Species_Name/Maker` in `maker_run.sh` shown below to your MAKER working directory.
```bash
#!/bin/bash
module load gffread
module load MakerP

while read i
do

cd /path/to/Species_Name/Maker/${i}

nohup maker maker_opts.ctl maker_bopts.ctl maker_exe.ctl &

done < chrs.list
```

Use `maker_run.sh` to submit MAKER runs in each directory setup by `maker_split_and_run.sh`.
```bash
bash maker_run.sh
```

### 5. Merge and Cleanup
#### GFF merge
Edit Species_Name before running.
```bash
global_gff="Species_Name.gff"; printf "##gff-version 3\n" >${global_gff} ;tail -n +2 */*maker.output/*_datastore/*/*/*/*gff | awk -F"\t" 'NF==9 && ($3=="gene" || $3 =="CDS" || $3 =="mRNA" || $3 =="exon" || $3 == "five_prime_UTR" || $3 == "three_prime_UTR" || $3 =="tRNA" )' >>${global_gff}
```

#### Make first protein database to guide filtering:
```bash
gffread -S -y Species_Name.prots -g Species_Name.genome.fa Species_Name.gff 
samtools faidx Species_Name.prots
``` 

Pick smallest protein.
```bash
sort -r -k2 -n Species_Name.prots.fai | cut -f -2 | grep est2genome
sort -k2 -n Species_Name.prots.fai | cut -f -2 | grep protein2genome | awk '$2 < 27 {print $1}' | sed 's/-mRNA-[0-9]*//' | sort | uniq > exclude_lt27.list`
```

#### Filter out the smallest proteins:
*Needs pip install gff3.*
```bash
printf "##gff-version 3\n" > Species_Name.ge27.gff ; python /path/to/gff_filter.py -e exclude_lt27.list -g Species_Name.gff >> Species_Name.ge27.gff
```

#### Make CDSs and UTRs unique:
```bash
python /path/to/gffPrepare/validate_gff.py --gff Species_Name.ge27.gff > Species_Name.ge27_uniq.gff
```

#### Rename:
Used to rename gff to the Patterson Lab/IWGC common naming convention. _Eleusine indica_ = EleIn. _Eleusine indica_ glyphosate-resistant = EleInR. _Species_name_ = SpeNa.
```bash
python /path/to/renameGff.py -g Species_Name.ge27_uniq.gff -t SpeNa > SpeNa.v2.gff
```

#### Sort GFF:
```bash
/path/to/gff3sort/gff3sort.pl --precise --chr_order natural SpeNa.v2.gff > SpeNa.v2.sorted.gff`
```

## Functional Annotation Usage
Unlike with structural annotation, the functional annotation pipeline is completely contained within a few custom scripts. Once configured, the below command is the only one to run.

### Functional annotation script:
Get your species' ID from NCBI.
```bash
/path/to/Functional_Annotation_v4.sh -G Species_Name.genome.fa -g SpeNa.v2.sorted.gff -s SpeNa -t 102 -i NCBI_Species_Name_ID -n Scientific_Name -c Common_Name
```
