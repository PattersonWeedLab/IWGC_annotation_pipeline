# IWGC_annotation_pipeline
Current annotation pipeline of the International Weed Genomics Consortium ([IWGC](https://www.weedgenomics.org/))  
Authors: Dr. Nathan D. Hall and Dr. Eric L. Patterson (P.I.)

----------------

# Summary of Functional Annotation

## 1. Isoform Selection
* Isoforms are selected using [AGAT](https://github.com/NBISweden/AGAT) which 
selects the longest isoforms with longest intact open read frame.
  
* Isoforms are then extracted from [GFF3](https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md)
  format using using [gffread](https://github.com/gpertea/gffread)
  
* These predicted proteins are used for all downstream analysis.

## 2. Location Prediction
* [MultiLoc2 Workstation Edition](https://github.com/NDHall/MultiLoc2-1/tree/workstation-edition)
  (ML2) is used predict protein location.
within the cell. Briefly, ML2 uses a trained machine learning classifer
  to predict protein location. For more on this program see the excellent 
  work by [Blum et al. (2009)](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-10-274)

## 3. InterPro Analysis
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
  
## 4. Homology Search
This search is for the best hit in each database using MMSeqs2. Which is an 
extremley fast method of searching that can be up to 36 times faster than BLAST [(Steinegger and Söding 2017)](https://www.nature.com/articles/nbt.3988).
### UniRef50
Is a clustered database that is regularly updated. See work by [Suzek et al. (2015)](https://academic.oup.com/bioinformatics/article/31/6/926/214968).
We used [UniRef_50 version 2021_03](https://ftp.uniprot.org/pub/databases/uniprot/previous_releases/release-2021_03/relnotes.txt)
and downloaded representative plants from the same release to add descriptions.
Cluster information and UniRef IDs were reported
## NCBI
All available proteins from fully sequenced plants were downloaded
May 2021 and put into a MMeqs2 database using [NCBI Datasets](https://www.ncbi.nlm.nih.gov/datasets/) and MMSeqs2 tool

## 5. Herbicide Interacting Genes
Herbicide interacting genes were identified from the literature and mapped to KEGG Ortholog. UniRef50 IDs
were converted to KEGG Orthologs and cross matched with herbicide interacting genes. Matches were flagged and 
Annotated with HRAC and WSSA classifiers.

## Output

The pipeline runs all analysis seperately, but a series of python scripts is used to collect, filter and 
and format the data. Results are reported in JSON format which is ideal for NOSQL databases and queries.
Given the number of Null values or missing values this format is ideal. The results are gappy because not every analysis 
returns a result for each isoform. This is expected. Finally, GFF3 is produced that contains Ontology Terms, and
Notes.

[fullsize poster](https://github.com/NDHall/NCWSS_2021/blob/main/media/NCWSS_poster_2021.pdf)

### Small Poster for Reference

![thumbnail](https://github.com/NDHall/NCWSS_2021/blob/main/media/NCWSS_poster.png)

----------------

# Installation

## Dependencies
The below program version numbers are the exact versions used by the current IWGC_annotation_pipeline. Different versions may also be compatible but no other versions have been verified to be compatible with this pipeline by the PattersonWeedLab.
* [RepeatModeler](https://github.com/Dfam-consortium/RepeatModeler) (version )
* [RepeatMasker](https://github.com/rmhubley/RepeatMasker) (version )
* [bedtools](https://github.com/arq5x/bedtools2) (version )
* [minimap2](https://github.com/lh3/minimap2) (version )
* [samtools](https://github.com/samtools/samtools) (version )
* [Cupcake](https://github.com/Magdoll/cDNA_Cupcake) (version )
* [gffread](https://github.com/gpertea/gffread/tree/master) (version )
* [maker](https://github.com/Yandell-Lab/maker) (version )
* [gff3](https://pypi.org/project/gff3/) (version 1.0.1)
* [BCBio](https://pypi.org/project/bcbio-gff/) (version 0.7.0)
* [gfftools](https://github.com/ihh/gfftools) (version )
* [AGAT](https://github.com/NBISweden/AGAT) (version )
* [interproscan](https://github.com/ebi-pf-team/interproscan) (version )
* [MMseqs2](https://github.com/soedinglab/MMseqs2) (version )
* [MultiLoc2 Workstation Edition](https://github.com/NDHall/MultiLoc2-1/tree/workstation-edition) (version )
