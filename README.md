# tractmapper

Converts a tract atlas to subject space and outputs the terminal coordinates of the fibers grouped by the tract.

## Requires:
wm_register_to_atlas_new.py -   https://github.com/SlicerDMRI/whitematteranalysis
                                module load whitematteranalysis/latest

tractconverter.py           -   https://github.com/MarcCote/tractconverter
                                module load tractconverter/latest

MIRTK.img                   -   https://github.com/TIGRLab/mirtk_singularity
                                singularity container

## Usage:

`get_subject_tract_coordinates.py --help`
