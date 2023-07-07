#!/bin/bash

print_usage() {
  printf "\


  Usage: Functional_prepare.sh -g <gff> -G <genome.fasta>
  -s <output base string> -t <threads>


  Variables:
  -G    Genome file path should be fasta this will be used for gffread step
  -g    gff path this will be used to get longest isoforms it will not be
        not modified with
  -s    this will be stem used for  all new files created. It should be
        informative, e.g. species name Ipomea_purpurea or code Ipopu
  -t    threads number of processors to use some processes will be threaded
        at the level of software while other processes, e.g. MultiLoc2 are
        threaded by virutue of file splitting starting specified number of
        processes in the background. This is done to produce a faster Multiloc2
        run.

  "
}

check_var () {
if [ -z ${1} ]; then print_usage ; exit 1 ; fi

}

while getopts 'g:G:k:t:s:' flag; do
  case "${flag}" in
    g) gff="${OPTARG}" ;;
    G) genome="${OPTARG}" ;;
    t) threads="${OPTARG}" ;;
    s) stem="${OPTARG}" ;;
    *) print_usage
       exit 1 ;;
  esac
done

check_var ${gff}
check_var ${genome}
check_var ${threads}
check_var ${stem}

printf "

Working dir: ${PWD}


Variables loaded and ready to go!

-g <gff>        ${gff}
-G <genome>     ${genome}
-s <stem>       ${stem}
-t <threads>    ${threads}


"




#-------------------------------------------------------------------------------

#   Config variables and functions are imported from pipeline_config.sh
#   pipeline_config.sh is sourced from the same folder as

#-------------------------------------------------------------------------------

echo $0 PATH to script

pipeline_config="$(echo ${0} | rev| cut -d"/" -f 2- |rev   )/pipeline_config.sh"

if [ -f ${pipeline_config} ] ; then
  source $( echo "$(echo ${0} | rev| cut -d"/" -f 2- |rev   )/pipeline_config.sh" )
else
  printf "\n\n${pipeline_config} is missing it is required for running this script\n\n"
  exit 1
fi



#-------------------------------------------------------------------------------

#   Create a single set of isoforms using agat

#-------------------------------------------------------------------------------

#-------------------------------output------------------------------------------
aa="${wd}/longest_dir/${stem}.aa"
cds="${wd}/longest_dir/${stem}.cds"

#   If no output file detected execute code block

if [ ! -f ${aa} ] ; then

    agat_start


    make_dir "longest_dir"

    module load ${agat_module}
    agat_sp_keep_longest_isoform.pl -gff ${gff} -o longest_dir/${stem}.gff

    module unload ${agat_module}
    #change to module load once gffread installed
    #add ${conda_module} to pipeline_config.sh
    load_conda ${conda_env001}

    gffread longest_dir/${stem}.gff -g ${genome} -J  -S -y ${aa}
    gffread longest_dir/${stem}.gff -g ${genome} -J  -S -x ${cds}

    conda deactivate


fi

#-------------------------------------------------------------------------------

#   Get MMSEQs hits in lieu of blast

#-------------------------------------------------------------------------------
module load ${mmseqs_module}
mmseqs_start
make_dir "mmseqs"

ncbi="mmseqs/${stem}_ncbi_select.ofmt6"

if [ ! -f ${ncbi} ] ; then
    mmseqs easy-search \
        ${aa} \
        ${ncbi_select_db} \
        ${ncbi} \
        "mmseqs/tmp" \
      --threads ${threads}  \
      --format-mode 0 \
      --start-sens 2 -s 7 \
      --sens-steps 3
    # remove  intermediate dir

    rm -rf "mmseqs/tmp"
fi

uniref="mmseqs/${stem}_uniref50.ofmt6"

if [ ! -f ${uniref} ]; then
    mmseqs easy-search \
        ${aa} \
        ${uniref_50} \
        ${uniref} \
        "mmseqs/tmp" \
      --threads ${threads}  \
      --format-mode 0 \
      --start-sens 2 -s 7 \
      --sens-steps 3

    # remove  intermediate dir

    rm -rf "mmseqs/tmp"
fi

herbres="mmseqs/${stem}_herbres.ofmt6"

if [ ! -f ${herbres} ]; then

    mmseqs easy-search \
        ${aa} \
        ${herbres_db} \
        ${herbres} \
        "mmseqs/tmp" \
      --cov-mode 0 \
       -c 0.30       \
      --threads ${threads}  \
      --format-mode 0 \
      --start-sens 2 -s 7 \
      --sens-steps 3

    # remove  intermediate dir

    rm -rf "mmseqs/tmp"
fi



#---------------------------Get Best Hit Leave Marker File----------------------

ncbi_best="mmseqs/${stem}_ncbi_best_hit.out"

if [ ! -f ${ncbi_best} ] ; then

    best_scoring_of6_blasthit ${ncbi} ${ncbi_best} && \
    rm ${ncbi} && \
    touch ${ncbi}

fi


uniref_best="mmseqs/${stem}_uniref_best_hit.out"

if [ ! -f ${uniref_best} ] ; then

    best_scoring_of6_blasthit ${uniref} ${uniref_best} && \
    rm ${uniref} &&\
    touch ${uniref}

fi

herb_hits="mmseqs/${stem}_filtered_herb_hits.tsv"
if [ ! -f ${herb_hits} ]; then

    awk ' $11 < 1e-5 {OFS="\t"; print $1,$2}' ${herbres} > ${herb_hits}

fi

module unload ${mmseqs_module}



#-------------------------------------------------------------------------------

#   Get Interpro results

#-------------------------------------------------------------------------------
ipr="${PWD}/iprscan/${stem}.iprscan.tsv"

if [ ! -f ${ipr} ] ; then

    iprscan5_start


    module load ${iprscan_module}

    make_dir iprscan
    interproscan.sh  \
              --formats GFF3 TSV \
              --goterms \
              --pathways \
              --iprlookup \
              --input ${aa} \
              --cpu ${threads} \
              --output-file-base iprscan/${stem}.iprscan

    module unload ${iprscan_module}
fi



#-------------------------------------------------------------------------------

#   Get Multiloc2 results

#-------------------------------------------------------------------------------

make_dir multiloc2 
module load singularity/3.9.7

go_list="${PWD}/multiloc2/${stem}_multilocGo.out"
prot_list="${PWD}/multiloc2/${stem}_prot.list"
simple_out="${PWD}/multiloc2/${stem}_ml2.txt"


multiloc2_start
if [ ! -f ${simple_out} ]; then
    splits="${PWD}/multiloc2/splits"
    make_dir ${splits}

    #---------------------------- Split --------------------------------------------

    awk -F"\t" '$14 ~/GO/ {print $1,$14}' ${ipr} | \
         sed  's/|/ /g' > ${go_list}

    egrep \> ${aa} | cut -b 2- | awk '{print $1}' > ${prot_list}

    cd ${splits} # this shortens output a little

    split -n l/${threads} ${prot_list}
    # split files will be labeled x[]a-z][a-z]

    declare -a pid_array
    for x in x[a-z][a-z]
    do
        mkdir ${x}_dir
        ml2_one_by_one.sh ${x} ${aa} ${go_list}  ${x}_dir &
        pid_array=( "${pid_array[@]}" "$!" )


    done
    pid_wait

    cat x*/ml2.out > ${simple_out}   
    cd - 


fi
