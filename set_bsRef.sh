#!/bin/bash -l

# Create context files for a reference genome

if [ "$#" != "1" ]; then
        echo "Usage: set_bsRef.sh [Ref]"
        exit 0
fi

#Load Modules
. /opt/shared/Modules/3.2.7/init/bash
module load bismark
module load bowtie/2-2.1.0

dir=$(mkdir -p ${basename $1})
cp $(readlink -f $1) $dir/

bismark_genome_preparation --bowtie2 $dir
