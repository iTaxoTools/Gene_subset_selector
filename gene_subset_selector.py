### gene_subset_selector.py####
### David Leisse ###
### david.leisse@uni-bielefeld.de ###

import copy
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

def get_mean_support(treefile: str, log_messages: list):
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
        log_messages.append(str(os.path.basename(treefile)) + " mean support: " + str(mean_support))
    else: 
        log_messages.append("!!!" + str(os.path.basename(treefile)) + " did not give any support values!")
        mean_support = 0.0

    return mean_support, log_messages


def top_n_files(filenumber: int, data: list, infolder: str, outfolder: str, log_messages: list):
    """
    Getting the top amount of files that have the highest scores.
    Copying and sorting Files from the input directory to the output directory by their corresponding files.

    :param Integer filenumber: Amount of top files that are selected
    :param List data: List of data containing filenames and corresponding scores
    :param String infolder: Path to input directory
    :param String outfolder: Path to output directory

    :return selected: List of selected top files 
    """

    tre_files = glob.glob(os.path.join(infolder, "Trees/*"))
    fasta_files = glob.glob(os.path.join(infolder, "Fasta/*"))
    if not fasta_files:
        sys.exit("Directory with .fasta files 'Fasta' not found or empty")
    
    all_file_list, missmatch_list = Utils.get_corresponding_files(tre_files, fasta_files = fasta_files)

    iter_list = copy.deepcopy(all_file_list)
    for basename, support in data:
        for idx,[file_a, file_b] in enumerate(iter_list):
            if basename in file_a or basename in file_b:
                all_file_list[idx].append(support)

    all_file_list.sort(reverse=True, key =lambda x: x[2])

    selected = all_file_list[:filenumber]
    not_selected = all_file_list[filenumber:]

    fasta_out_selected = os.path.join(outfolder,"Selected_alignments")
    tre_out_selected = os.path.join(outfolder,"Selected_trees")
    fasta_out_not_selected = os.path.join(outfolder,"Not_selected_alignments")
    tre_out_not_selected = os.path.join(outfolder,"Not_selected_trees")

    for entry in all_file_list:
        if entry in selected:
            shutil.copy(entry[0], tre_out_selected)
            shutil.copy(entry[1], fasta_out_selected)

        elif entry in not_selected:
            shutil.copy(entry[0], tre_out_not_selected)
            shutil.copy(entry[1], fasta_out_not_selected)

    for entry in missmatch_list:
        log_messages.append(str(entry)+ " no corresponding file!!")

    return selected, log_messages

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

    log_messages = []

    #### Getting scores #######
    scores = OrderedDict()
    if criterion == "support":
        log_messages.append("Criterion selected: Tree support")
        tre_files = glob.glob(os.path.join(input, "Trees/*"))
        if not tre_files:
            sys.exit("!! Directory with .tre files 'Trees' not found or empty")
        for tre in tre_files:
            mean_support, log_messages = get_mean_support(tre, log_messages)
            name = os.path.splitext(os.path.basename(tre))[0]
            scores[name] = mean_support

    elif criterion == "alignment":
        log_messages.append("Criterion selected: Alignment")
        print("in work")

    elif criterion == "missingdata":
        log_messages.append("Criterion selected: Missingdata")
        print("in work")

    elif criterion == "ntaxa":
        log_messages.append("Criterion selected: N Taxa")
        print("in work")

    #### Selecting #########
    sorted_scores = sorted(scores.items(), key=lambda x: x[1], reverse=True)
    print(sorted_scores)
    if filenumber:
        log_messages.append("Selection chosen: filenumber")
        selected_scores, log_messages = top_n_files(filenumber, sorted_scores, input, output, log_messages)
        coloumns = ["file",criterion + " score"]
        sorted_scores.insert(len(selected_scores),["\nNot selected files:",""])
        Utils.write_tsv_output(os.path.join(output, "file_selection.tsv"), coloumns, sorted_scores)

    elif percentage:
        log_messages.append("Selection chosen: percentage")
        print("in work")

    elif cutoff:
        log_messages.append("Selection chosen: cutoff")
        print("in work")
    
    #### Log file 
    Utils.write_log_file(log_messages, output)

__Main__("--dir /Users/david/Development/Gene_subset_selector-1/test_subset --out /Users/david/Development/Gene_subset_selector-1/test_subset --files 3 --crit support")


"""
if "--dir" in sys.argv and "--out" in sys.argv and "--crit" in sys.argv:
    __Main__(sys.argv)

else:
    print(__usage__)
    """