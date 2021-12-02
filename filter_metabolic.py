#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import argparse
import os

""" compared to main branch:
    - takes abundances file without the empty row above.
    - exports sheets as .tsv rather than .csv
"""

parser = argparse.ArgumentParser()

parser.add_argument("abundances", help = "file with MAG abundance per sample") # file for abundances
parser.add_argument("-f", "--folder", default = os.getcwd(), help = "folder to look for METABOLIC output worksheets1-6. Defaults to CWD") 
parser.add_argument("-t", "--threshold", default = 1.0, type=float, help = "Minimum abundance in at least 1 sample. Default: 1.0")
parser.add_argument("-o", "--output", default ="", help = "specify prefix for destination folder that will be created") # set a prefix for the output folder
parser.add_argument("-x", "--taxonomy", action="store_true", help = "if set, MAG name is replaced with associated taxonomy from the 'abundances' file")
args = parser.parse_args()



def fix_underscores(mags_column):
    return [name.replace(".", "_") for name in mags_column.to_list()] # replace all "." in the column with "_"

def foldername(name):
    copy = 0
    if name != "":
        name = name + "_"
    folder = name + "filter_METABOLIC"
    if os.path.exists(folder):
        copy = 1
        folder = name + f"filter_METABOLIC({str(copy)})"
        while os.path.exists(folder):
            copy += 1
            folder = name + f"filter_METABOLIC({str(copy)})"
    return folder
    
def log():
    message = ["METABOLIC folder:\t" + os.path.realpath(args.folder) + "\n",
            "MAGs abundance file:\t" + os.path.realpath(args.abundances) + "\n",
               "Samples:\t" + ", ".join(abundances_cols.columns.tolist()) + "\n",
               "Threshold:\t" + str(args.threshold) + "\n",
               "Output folder:\t" + os.path.realpath(output_folder)]
    return message

def create_log():
    with open(os.path.join(output_folder,"log.txt"), "w") as f:
        for line in log():
            f.write(line)        

def process_worksheet(sheet):
    """per quanto riguarda la necessità di tenere solo un pattern di colonne,
    come nel primo foglio solo le colone hit numbers, invece di fargli
    fare lo slicing su mag_hits[None::None] in ogni caso, ha più senso
    dargli un parametro "Needs_slicing" == True per eseguire lo slicing
    solo se serve"""
    
    if sheet == 1:
        to_strip = [" Hmm presence", " Hit numbers", " Hits"]
        mags_starting_column = 10
        pattern = (1, 3)
        keeps_present = False
        keeps_nonzero = True
        sheet_name = "HMMHitNum.tsv"
    if sheet == 2:
        to_strip = [" Function presence"]
        mags_starting_column = 3
        pattern = (None, None)
        keeps_present = True
        keeps_nonzero = False
        sheet_name = "FunctionHit.tsv"
    if sheet == 3:
        to_strip = [" Module presence"]
        mags_starting_column = 3
        pattern = (None, None)
        keeps_present = True
        keeps_nonzero = False
        sheet_name = "KEGGModuleHit.tsv"
    if sheet == 4:
        to_strip = [" Module step presence"]
        mags_starting_column = 4
        pattern = (None, None)
        keeps_present = True
        keeps_nonzero = False
        sheet_name = "KEGGModuleStepHit.tsv"
    if sheet == 5:
        to_strip = [" Hit numbers", " Hits"]
        mags_starting_column = 1
        pattern = (0, 2)
        keeps_present = False
        keeps_nonzero = True
        sheet_name = "dbCAN2Hit.tsv"
    if sheet == 6:
        to_strip = [" Hit numbers", " Hits"]
        mags_starting_column = 1
        pattern = (0, 2)
        keeps_present = False
        keeps_nonzero = True
        sheet_name = "MEROPSHit.tsv"

        
    try:
        work = pd.read_csv(os.path.join(args.folder, f"METABOLIC_result_worksheet{sheet}.tsv"), sep ="\t")
    except FileNotFoundError:
        print("**************ERROR: could not find the worksheet: " + os.path.join(args.folder, f"METABOLIC_result_worksheet{sheet}.tsv"), sep ="\t")
    else:
        info = work.drop(work.columns[mags_starting_column:], axis=1) # dataframe with reaction columns only
        hits = work.drop(work.columns[0:mags_starting_column], axis=1) # dataframe with MAG hits only
        
        for string in to_strip:
            hits.columns = hits.columns.str.replace(string, "")
        
        colonne = hits.columns.isin(to_search) # boolean array to select only the columns of the MAGs specified in "to_search"
        mag_hits = hits.loc[:, colonne]     # dataframe with the aforementioned columns
        
        #mag_hits.to_csv("by_worksheet_new_work1_hits.csv", sep="\t", header=True, index = False)
        mag_hits = mag_hits.iloc[:, pattern[0]::pattern[1]] # for each MAG, only keep the appropriate column
           
        if args.taxonomy == True and missing_taxonomy == False:
            tax_hits = mag_hits.rename(axis="columns", mapper=taxonomy) # replace MAGs name with their taxonomy
        else:
           tax_hits = mag_hits
        
        new_sheet = pd.concat([info, tax_hits], axis = 1) # join the dataframe of the reactions of the selected MAGs.
        #new_sheet.to_csv("new_work1.csv", sep="\t", header=True, index = False)
           
        if keeps_present:
            only_present = new_sheet[(new_sheet.iloc[:, mags_starting_column::] == "Present").any(axis=1)]  # only the rows with at least 1 "Present" in the MAGs hits columns  
        elif keeps_nonzero:
            only_present = new_sheet[(new_sheet.iloc[:, mags_starting_column::] != 0).any(axis=1)]   # only the rows with at least 1 non-zero value in the MAG hits columns
        
        output_name = f"filter_METABOLIC_{sheet_name}"
        only_present.to_csv(os.path.join(output_folder, output_name), sep = "\t", header = True, index = False)

############################### MAIN
if os.path.isfile(args.abundances) == False:
    print("*********Could not find the specified abundances file: " + str(args.abundances))
    raise SystemExit


df = pd.read_csv(args.abundances, sep='\t', header=0)
df["user_genome"] = fix_underscores(df["user_genome"])

if args.taxonomy == True:
    try:
        taxonomy = dict(df[["user_genome","taxonomy"]].values) # a dictionary with MAG id and taxonomy from the provided abundances file.
    except KeyError:
        # if the file has no "taxonomy" column, print a warning and move on, skip the taxonomy replacement
        print("\n *********Could not replace MAG names with taxonomy: your abundances csv file does not contain the 'taxonomy' column")
        missing_taxonomy = True
    else:
        missing_taxonomy = False
        for magname, tax in taxonomy.items():
            # if a mag has no assigned taxonomy, leave the MAG id
            if pd.isna(tax):
                taxonomy[magname] = magname
    
        
output_folder = foldername(args.output)
os.mkdir(output_folder)

abundances_cols = df.iloc[:, 2:] # we look for MAG abundances from third column onwards
small_df = df[(abundances_cols >= args.threshold).any(axis = 1)]
#small_df.to_csv(os.path.join(args.output, f"plus_{args.output}_tax_select_abundant_mags.csv"), sep='\t', header=True, index=False)
small_df.to_csv(os.path.join(output_folder, "abundant_mags.csv"), sep='\t', header=True, index=False)


to_search = small_df["user_genome"].tolist()
to_search.sort() # because MAGs are sorted alphabetically in the METABOLIC worksheet

for i in range(1,7):
    process_worksheet(i)
    
create_log()
#print("\n")
print(*log())
#print("\n")
