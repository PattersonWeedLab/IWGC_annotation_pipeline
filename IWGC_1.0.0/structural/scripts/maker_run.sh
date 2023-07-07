#!/bin/bash

conda activate iwgcsa001

module load gffread

module load MakerP

while read i
do

cd /data/projects/01_struc_anno/results/Eluesine_indica_Gly_R/maker/${i}

nohup maker maker_opts.ctl maker_bopts.ctl maker_exe.ctl &

done < chrs.list
