# Current annotation pipeline of the International Weed Genomics Consortium ([IWGC](https://www.weedgenomics.org/))  
**Authors: [Dr. Nathan D. Hall](https://github.com/NDHall) (Developer), [Nicholas A. Johnson](https://github.com/Scrumpis) (Documentation), and Dr. Eric L. Patterson (P.I.)**  
**GitHub Repo: [IWGC_annotation_pipeline](https://github.com/PattersonWeedLab/IWGC_annotation_pipeline)**  



# Introduction
  
## Summary of Structural Annotation

### 1. Handling Repeats
* Repeat regions are annotated using RepeatModeler (links and versions of all programs in [Dependencies](#Dependencies) section).
  
* Annotated repeat regions are masked with RepeatMasker to reduce
  the computational burden of further analysis.

* bedtools is used to soft mask the genome using the output
  of RepeatMasker.
  
### 2. Map Isoseq Reads to Masked Genome
#### minimap2 alignment
* IsoSeq reads are mapped to the repeat-masked genome using minimap2.

#### SAMtools SAM to BAM
* SAMtools is used to convert the minimap2 alignment SAM file into a BAM file.

### 3. Collapse Isoforms
#### CupCake:
* Isoforms are collapsed with CupCake.  

### 4. Maker
* Genome, collapsed cDNA from Cupcake, repeat libraries from RepeatModeler, and a protein FASTA from a close relative species are fed into MAKER

### 5. Merge and Cleanup
* Genes that produced proteins under 32 amino acids long were removed from further annotation with only the longest proteins from each gene and unique untranslated regions (UTRs) used for functional annotation.


## Summary of Functional Annotation

### 1. Isoform Selection
* Isoforms are selected using [AGAT](https://github.com/NBISweden/AGAT) which 
selects the longest isoforms with longest intact open read frame.
  
* Isoforms are then extracted from [GFF3](https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md)
  format using using [gffread](https://github.com/gpertea/gffread)
  
* These predicted proteins are used for all downstream analysis.

### 2. Location Prediction
* [MultiLoc2 Workstation Edition](https://github.com/NDHall/MultiLoc2-1/tree/workstation-edition)
  (ML2) is used predict protein location.
within the cell. Briefly, ML2 uses a trained machine learning classifer
  to predict protein location. For more on this program see the excellent 
  work by [Blum et al. (2009)](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-10-274)

### 3. InterPro Analysis
* InterPro is run locally using iprscan5 version 5.47-82.0
* It uses several programs to predict InterPro Accessions which 
can be broken down into the following categories
    
    1. Sites
        * Conserved Site
        * Active Site
        * Binding Site
        * PTM
    2. Repeats
    3. Domains
    4. Families
    5. Homologous SuperFamilies
       
* For more information see release notes [here](https://www.ebi.ac.uk/interpro/release_notes/84.0/)
* IPRSCAN Results provide [MetaCyc](https://metacyc.org/) accessions but not descriptions. To obtain descriptions
a custom python script was used to extract Pathway IDs and link them to their 
  descriptions. MetaCyc is a proprietary database that is made conditionally available to academic researchers 
  free of charge.
  
### 4. Homology Search
This search is for the best hit in each database using MMSeqs2. Which is an 
extremley fast method of searching that can be up to 36 times faster than BLAST [(Steinegger and Söding 2017)](https://www.nature.com/articles/nbt.3988).
#### UniRef50
Is a clustered database that is regularly updated. See work by [Suzek et al. (2015)](https://academic.oup.com/bioinformatics/article/31/6/926/214968).
We used [UniRef_50 version 2021_03](https://ftp.uniprot.org/pub/databases/uniprot/previous_releases/release-2021_03/relnotes.txt)
and downloaded representative plants from the same release to add descriptions.
Cluster information and UniRef IDs were reported
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
BuildDatabase -name species_name genome.fasta
``` 

Next use repeat database from last step to model repeats. -pa 1 uses 4 threads, -pa 25 uses 100 threads. Adjust accordingly. 
```bash
RepeatModeler -database species_name -pa 25 -LTRStruct
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
samtools view -b -T path/to/og_genome.fasta minimap2_alignment.sam > minimap2_alignment.bam
```

### 3. Collapse Isoforms with Cupcake
#### CupCake:
Collapse isoforms with CupCake.
```bash
python path/to/collapse_isoforms_by_sam.py --input ISOseq.fq --fq -b minimap2.sorted.bam -o output_genome_rootname
```

### 4. Maker


### 5. Merge and Cleanup
####
```bash
global_gff="EleIndGlyRes01.gff"; printf "##gff-version 3\n" >${global_gff} ;tail -n +2   Chr*/*maker.output/*_datastore/*/*/Chr*/*gff | awk -F"\t" 'NF==9 && ($3=="gene" || $3 =="CDS" || $3 =="mRNA" || $3 =="exon" || $3 == "five_prime_UTR" || $3 == "three_prime_UTR" || $3 =="tRNA" )' >>${global_gff}
```

#### Make first protein database to guide filtering:
```bash
gffread -S -y EleIndGlyRes03.prots -g Eindica_glyres.genome.fa EleIndGlyRes03.gff 
samtools faidx EleIndGlyRes03.prots
``` 

Pick smallest protein.
```bash
sort -r -k2 -n genome_name.prots.fai | cut -f -2 | grep est2genome
sort -k2 -n genome_name.prots.fai | cut -f -2 | grep protein2genome | awk '$2 < 27 {print $1}' | sed 's/-mRNA-[0-9]*//' | sort | uniq > exclude_lt27.list`
```

#### Filter them out:
*Needs pip install gff3.*
```bash
printf "##gff-version 3\n" > EleIndGlyRes03.ge27.gff ; python /data/projects/01_struc_anno/workflow/gff_filter.py -e exclude_lt27.list -g EleIndGlyRes03.gff >> EleIndGlyRes03.ge27.gff
```

#### Make CDSs and UTRs unique:
```bash
python /data/projects/iwgc_annotation/workflow/keyGene/src/Key_Gene_Scripts/gffPrepare/validate_gff.py --gff EleIndGlyRes03.ge27.gff > EleIndGlyRes03.ge27_uniq.gff
```

#### Rename:
```bash
python /data/projects/01_struc_anno/workflow/renameGff.py -g EleIndGlyRes03.ge27_uniq.gff -t EleInR > EleIndR02.gff
```

#### Sort GFF:
```bash
/path/to/gff3sort/gff3sort.pl --precise --chr_order natural EleIndR02.gff > EleIndR02.sorted.gff`
```

## Functional Annotation Usage
Unlike with structural annotation, the functional annotation pipeline is completely contained within a few custom scripts. Once configured, the below command is the only one to run.

### Functional annotation script:
```bash
/path/to/Functional_Annotation_v4.sh -G genome_name.genome.fa -g genome_name.sorted.gff -s EleInR -t 102 -i -n scientific_name -c common_name
```
