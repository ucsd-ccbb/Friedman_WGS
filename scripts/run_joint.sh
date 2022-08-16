#!/bin/bash

metadata=/shared/workspace/projects/friedman/metadata/friedman_joint_metadata_04182022.csv
workspace=/scratch/joint

qsub \
-v workspace=$workspace \
-v metadata=$metadata \
/shared/workspace/projects/friedman/scripts/joint_genotyping.sh
