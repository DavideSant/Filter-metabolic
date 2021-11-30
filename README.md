# Filter-metabolic
This script filters the output of METABOLIC to extract features present in abundant MAGs only.
1. Select MAGs with an abundance higher than the user-specified threshold.
2. Select the features present in at least 1 of the abundant MAGs.

### Requirements
- pandas library for python. OS, argparse modules from the standard library.

### Input
- The abundances file.
- The 6 separate .tsv worksheets, output of METABOLIC.

#### Abundances file
_Please use the provided template file and fill it in with your data._
A **tab-separated** .csv file with the abundance of each MAG in the samples.
- You can name the sample columns with the name of your samples. 
- A MAG will be selected if it's more abundant than threshold in at least 1 sample. 
- You can remove unused sample columns, and you can add other sample columns if you have more than 2 samples.
- In the taxonomy column you can provide taxonomy associated to each MAG.

### Output
The files are created into the folder `<prefix>_filter_METABOLIC`.
- `abundant_mags.csv`: .csv file with the abundant MAG only.
- The 6 filtered .tsv worksheets. Optionally, the MAG id will be replaced with the associated taxonomy provided in the abundances file.

### Usage
`$ ./filter_metabolic.py <MAG abundances .csv> -f <folder with METABOLIC worksheets> -t <abundance threshold> -o <output folder prefix> -x`
- `<MAG abundances .csv>`: is the only required argument. It's the file with the abundances of MAGs.
- `<folder with METABOLIC worksheets>`: is the folder where the METABOLIC worksheets are stored. Defaults to CWD.
- `<abundance threshold>`: the minimum abundance for a MAG to be selected as abundant. Defaults to 1.
- `<output folder prefix>`: prefix added to "filter_METABOLIC" name for the folder containing the output.
- `-x`: if set the script will perform the taxonomy replacement.
