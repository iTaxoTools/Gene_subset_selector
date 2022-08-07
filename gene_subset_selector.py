### gene_subset_selector.py####
### David Leisse ###
### david.leisse@uni-bielefeld.de ###

import Utils
import dendropy
from collections import OrderedDict
import os.path
import sys
import glob
from genericpath import isdir
from os import mkdir
import shutil

__usage__ = """

            python3 gene_subset_selector.py
            --dir <PATH_TO_INPUT_DIRECTORY>
            --out <PATH_TO_OUTPUT_DIRECTORY>
            --crit <CRITERION:support,alignment,missingdata,ntaxa>

            categories (select one):
            --files <INTEGER_AMOUNT_OF_TOP_FILES>
            --perc <INTEGER_PERCENTAGE_OF_TOP_FILES>
            --cutoff <FLOAT_SCORE_CUTOFF>

"""

def get_mean_support(treefile: str):
    """
    Calculating the mean support of nodes in a tree. 
    Either of newick or nexus formatted trees.

    :param String treefile: Path to treefile
    :return mean_support: Mean support value of corresponding tree
    """

    if os.path.splitext(treefile)[1] == ".tre":
        tree = Utils.load_tree(treefile,"path", "newick")
    
    elif os.path.splitext(treefile)[1] == ".nex":
        tree = Utils.load_tree(treefile,"path", "nexus")

    sum_support = 0.0
    n_nodes = 0
    for i in tree.__iter__():
        if i.label:
            sum_support += float(i.label)
            n_nodes += 1
    
    if n_nodes and sum_support:
        mean_support = sum_support/n_nodes
    else: 
        print(str(os.path.basename(treefile)) + " did not give any support values!")
        mean_support = 0.0

    return mean_support


def top_n_files(filenumber: int, data: list, infolder: str, outfolder: str):
    """
    Getting the top amount of files that have the highest scores.
    Copying and sorting Files from the input directory to the output directory by their corresponding files.

    :param Integer filenumber: Amount of top files that are selected
    :param List data: List of data containing filenames and corresponding scores
    :param String infolder: Path to input directory
    :param String outfolder: Path to output directory

    :return selected: List of selected top files 
    """

    selected = data[:filenumber]
    not_selected = data[filenumber:]

    tre_files = glob.glob(os.path.join(infolder, "Trees/*"))
    fasta_files = glob.glob(os.path.join(infolder, "Fasta/*"))

    for file, support in selected:
        fasta_out = os.path.join(outfolder,"selected_alignments")
        tre_out = os.path.join(outfolder,"selected_trees")

        fasta_file = [x for x in fasta_files if file in x][0]
        tre_file = [x for x in tre_files if file in x][0]

        fastafilepath = os.path.join(infolder,"Fasta",fasta_file)
        trefilepath = os.path.join(infolder,"Trees",tre_file)
        shutil.copy(fastafilepath, fasta_out)
        shutil.copy(trefilepath, tre_out)

    for file, support in not_selected:
        fasta_out = os.path.join(outfolder,"not_selected_alignments")
        tre_out = os.path.join(outfolder,"not_selected_trees")

        fasta_file = [x for x in fasta_files if file in x][0]
        tre_file = [x for x in tre_files if file in x][0]

        fastafilepath = os.path.join(infolder,"Fasta",fasta_file)
        trefilepath = os.path.join(infolder,"Trees",tre_file)
        shutil.copy(fastafilepath, fasta_out)
        shutil.copy(trefilepath, tre_out)

    return selected

def __Main__(args):

    if type(args) == str:
        args = args.split(" ")

    input = args[args.index("--dir")+1 ]
    output = args[args.index("--out") +1 ]
    criterion = args[args.index("--crit")+1 ]
    
    filenumber = 0
    percentage = 0
    cutoff = 0
    if "--files" in args:
        filenumber = int(args[args.index("--files")+1])
    elif "--perc" in args:
        percentage = int(args[args.index("--perc")+1])/100
    elif "--cutoff" in args:
        cutoff = float(args[args.index("--cutoff")+1])

    category = "nucleotide"
    if "--type" in args:
        category = args[args.index("--type")+1 ]

    outfolders = ["Selected_alignments","Selected_trees","Not_selected_trees","Not_selected_alignments"]
    for outfolder in outfolders:
        path = os.path.join(output, outfolder)
        if isdir(path):
            shutil.rmtree(path)
        mkdir(path)

    #### Getting scores #######
    scores = OrderedDict()
    if criterion == "support":
        tre_files = glob.glob(os.path.join(input, "Trees/*"))
        for tre in tre_files:
            mean_support = get_mean_support(tre)
            name = os.path.splitext(os.path.basename(tre))[0]
            scores[name] = mean_support

    elif criterion == "alignment":
        print("in work")

    elif criterion == "missingdata":
        print("in work")

    elif criterion == "ntaxa":
        print("in work")

    #### Selecting #########
    sorted_scores = sorted(scores.items(), key=lambda x: x[1], reverse=True)
    print(sorted_scores)
    if filenumber:
        selected_scores = top_n_files(filenumber, sorted_scores, input, output)
        coloumns = ["file",criterion + " score"]
        sorted_scores.insert(len(selected_scores),["\nNot selected files:",""])
        Utils.write_tsv_output(os.path.join(output, "file_selection.tsv"), coloumns, sorted_scores)

    elif percentage:
        print("in work")

    elif cutoff:
        print("in work")



if "--dir" in sys.argv and "--out" in sys.argv and "--crit" in sys.argv:
    __Main__(sys.argv)

else:
    print(__usage__)
