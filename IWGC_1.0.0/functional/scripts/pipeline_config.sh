#!/bin/bash

#--------------	Part I Variables --------------------
project_scripts="/data/projects/iwgc_annotation/workflow/scripts/"
agat_module="agat"
mmseqs_module="MMseqs2"
iprscan_module="iprscan/5.47-82.0"
singularity_module="singularity" #Multiloc2 is being run from a container
conda_env001="iwgcfa001"
ncbi_select_db="/data/projects/iwgc_annotation/resources/ncbi/seqs/current/GenomesDatabase.v1.prot.mmseqs"
uniref_50="/data/projects/iwgc_annotation/resources/uniref/mmseqs_db/uniref50"
herbres_db="/data/projects/iwgc_annotation/resources/herb_res_prots/total_KO_labeled.mmseqs_dbs"
#Path to Custom Multiloc2 Package
ml="/data/projects/iwgc_annotation/workflow/singularity-images/multiloc2_v3.img"
wd=${PWD}
ml2_path="${project_scripts}:/usr/share/conda/anaconda3/bin/:/usr/share/conda/anaconda3/bin/python3:/opt/software/MultiLoc2/1.0:/opt/software/MultiLoc2/1.0/src:/opt/software/signalp/5.0b/bin:/opt/software/signalp/5.0b:/opt/software/TMHMM/2.0c/bin:/opt/software/phobius/1.01:/opt/software/GNU/5.5.0/bin:/opt/software/GNU/5.5.0:/etc/alternatives/jre_11_openjdk/bin:/opt/software/iprscan/5.47-82.0:/opt/software/libsvm/3.25:/opt/software/BLAST+/2.3.0/bin:/opt/software/BLAST/2.2.26/bin:/usr/share/conda/anaconda3/bin:/usr/share/conda/anaconda3/condabin:/usr/share/Modules/bin:/usr/local/bin:/usr/bin:/usr/local/sbin:/usr/sbin:/usr/share/conda/anaconda3/condabin:/data/projects/iwgc_annotation/workflow/keyGene/shellScripts:/home/hallna12/.local/bin:/home/hallna12/bin"
ml2_python_path=""
#-------------- Part II Variables --------------------

pipeline_root="/data/projects/functional_annotation/"
package_root="/data/projects/iwgc_annotation/workflow/keyGene"



export PATH="${PATH}:${project_scripts}"
export PYTHONPATH=$PYTHONPATH:${package_root} # so python can find root modules

#-------------------------------------------------------------------------------

#   functions

#-------------------------------------------------------------------------------

make_dir () {
if [ ! -d ${1} ] ; then mkdir ${1} ; fi
}

load_conda () {
source ~/.bashrc
conda activate ${1}
}

pid_wait () {
for pid in ${pid_array}
do
wait ${pid}
done
pid_array=()

}

best_scoring_of6_blasthit () {
# sort and awk need to use tabs or it will sort based on the wrong column.
# The ncbi_select_v1 database has headers with spaces in the names and these are
# reported back.
sort -k12 -t $'\t' -nr ${1} | awk -F "\t" ' ! a[$1]++ && $11 < 1e-5' > $2
}

agat_start () {


printf "

Extract longest isoforms using agat and gffread


Input:
gff file    ${gff}
genome      ${genome}

Output:

outputdir      ${PWD}/longest_dir/

gff file       ${stem}.gff
seq file       ${stem}.aa

"


}

mmseqs_start () {
printf "
input:
seq file       ${stem}.aa

Output:

outputdir      ${PWD}/mmseqs/

ncbi select    ${stem}_ncbi_select.ofmt6
uniref_50      ${stem}_uniref50.ofmt6


"
}

iprscan5_start () {
printf "

Starting iprscan5 v5.47.82.0.Py3

Input:
seq file    ${aa}

Output:
ipr.tsv



"
}

multiloc2_start () {
printf "
Input:
go terms extracted from     ${ipr}
seq file                    ${aa}

Output:
go_list     ${go_list}
prot_list   ${prot_list}
simple_out  ${simple_out}

"
}
##====================================================================================================================


#	Part II (Mergeing Results)


##====================================================================================================================


