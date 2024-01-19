#!/bin/bash 
set -e
set -o pipefail

wd=$PWD

refwY=$1
refwoY=$2

PIPELINE=/data/Phillippy2/projects/rpc/99.codes
for sexChr in XY XX
do
if [ "$sexChr" == "XX" ] ; then
        ref=$refwoY
elif [ "$sexChr" == "XY" ] ; then
        ref=$refwY
fi
echo -e "indexing $ref"
if [[ ! -f $ref.bwt ]]; then
        cpus=4
        mem=10g
        name=index.$sexChr
        script=$PIPELINE/bwa/bwa_index.sh
        args="$ref"
        partition=quick
        walltime=4:00:00
        path=$wd
        log=$wd/logs/$name.%A.log

        mkdir -p $wd/logs

        echo "\
        sbatch -J $name --mem=$mem --partition=$partition --cpus-per-task=$cpus -D $path $extra --time=$walltime --error=$log --output=$log $script $args"
        sbatch -J $name --mem=$mem --partition=$partition --cpus-per-task=$cpus -D $path $extra --time=$walltime --error=$log --output=$log $script $args | awk '{print $NF}' > $wd/$sexChr.index.jid
        echo $index
fi
done

