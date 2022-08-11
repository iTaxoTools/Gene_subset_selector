# Gene_subset_selector
!!! Only **support** criterion and selection by **filenumber** has been implemented so far!!!

Selecting subset(s) of genes through different methods.
Alignment files can be compared and ranked through different methods. 
A certain amount of files is selected as the wanted subset of genes.

## Preparing data

One input directory is necessary. It should contain a folder called "Fasta" containing all necessary alignment files
in a fasta format and a folder called "Trees" that contains all corresponding tree files.
Each corresponding Tree and alignment files have to be named the same (except for the suffix).

## Running the script

```python3 gene_subset_selector.py --dir <Path to input directory> --out <Path to output directory> --crit <Selected criterion> --files <Amount of top files>```

`--dir` Defines the directory where all alignment and tree files can be found

`--out` Defines the directory where all output folders and files will be created

`--crit` Defines the criterion by which a score for each alignment file is calculated: support, alignment, missingdata, ntaxa

Select one of these selection methods in addition:

`--files` Defines the top amount of files that are supposed to be selected in the ranking.

`--perc` Defines the top percentage of files that are supposed to be selected in the ranking.

`--cutoff` Defines a cutoff for the score under which a file is no longer selected.

## Criterions

**support** Calculates a support value score by going through each corresponding tree file and calculating its mean support value.

**alignment**

**missingdata**

**ntaxa**

## Output
The script creates 4 folders in the output directory:

**Selected alignments** contains all fasta alignment files that have been selected by the chosen criterion.

**Selected trees** contains all tre files corresponding to the selected fasta alignment files.

**Not selected alignments** contains all fasta alignment files that were not sufficiant to the chosen criterion.

**Not selected trees** contains all tre files corresponding to the not selected fasta alignment files.
