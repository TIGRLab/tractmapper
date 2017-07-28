#!/bin/bash
atlasFile=/scratch/twright/tractSOP/step_2_output/clustered_whole_brain.vtp
clusterDir=/scratch/twright/tractSOP/step_2_output/
mrmlFile=/scratch/twright/tractSOP/step_2_output/clustered_tracts_display_100_percent_aem.mrml
subjectFile=/archive/data/SPINS/pipelines/dtiprep/SPN01_CMH_0004_01/SPN01_CMH_0004_01_02_DTI60-1000_05_Ax-DTI-60plus5_SlicerTractography.vtk
anatFile=/archive/data/SPINS/data/nii/SPN01_CMH_0004_01/SPN01_CMH_0004_01_02_DTI60-1000_05_Ax-DTI-60plus5.nii.gz
workdir=/tmp/work/
container=/archive/code/containers/MIRTK/MIRTK.img

./get_subject_tract_coordinates.py --mirtk_file=$container $atlasFile $clusterDir $mrmlFile $subjectFile $anatFile
