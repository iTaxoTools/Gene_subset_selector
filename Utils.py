##### Utils.py #####
##### David Leiße #####
##### david.leisse@uni-bielefeld.de #####

import dendropy
from os import mkdir
from genericpath import isdir
import os.path
import itertools
import re

def load_tree(string: str, mode: str, schema: str) -> dendropy.Tree:
    """
    Loading tree from .tre file and plotting it.
    Working with Newick schema until now.

    :param String string: Relative or absolute path to .tre file or data string
    :param String mode: Choosing between 'path' and 'string'
    :param String schema: I.e. 'newick' or 'nexus' depending on .tre file

    :return tree: Dendropy.Tree object
    """
    if mode == "path":
        tree = dendropy.Tree.get(path = string, schema = schema,rooting ="force-unrooted")
    elif mode == "string":
        tree = dendropy.Tree.get(data = string, schema = schema)
        
    return tree

def load_fasta_ali_file(file: str) -> list[str,str]:
    """
    Loading fasta file into list. List containing lists with each an ID and sequence

    :param String file: full path to file
    :return ls: List containing lists with each an ID and sequence  
    """

    with open(file, "r") as f:
        line = f.readline()
        characters = ["@"," "]
        sequence = []
        header = ""
        ls = []
        while line:
            
            while line and not line.startswith(">"):
                if not sequence == "\n":
                    sequence.append(line.strip())
                line = f.readline()

            if sequence and header:
                ls.append([header, "".join(sequence)])
                
            sequence = []
            header = line[1:].strip()
            for character in characters:
                if character in header:
                    header.replace(character, "_")

            line = f.readline()
    
    return ls

def write_decont_output(directory: str, file: str, seq_ls: list[str,str], type: str):
    """
    Writing output files from modified dictionary

    :param List mod_ali_ls: modified list with contents for decontaminated .ali files
    :param String ali_file: path to original .ali file
    """

    output_directory = directory + "/decontaminated"
    if not isdir(output_directory):
        mkdir(output_directory)

    file = output_directory + "/" + file

    with open(file, "w") as out:
        if type == "protein":
            out.write("#\n#\n")
        for id,seq in seq_ls:
            out.write(">" + id + "\n")
            out.write(seq + "\n")

def get_index_in_list(value: any, list: list) -> int:
    """
    Getting index of list of lists by value.

    :param Any value: value that is searched in list of lists
    :param List list: corresponding List of lists

    :return index: Index that corresponds to list that contains value
    """
    index = -1
    for idx,x in enumerate(list):
        if value in x:
            index = idx

    if index == -1:
        return None

    else:
        return index

def get_corresponding_files(tre_files: list, ali_files:list = [], fasta_files: list = []) -> list:
    """
    Getting all corresponding files to a .tre file. 
    Corresponding files are name similar to .tre file.
    Usually .fasta / .ali files

    :param List tre_files: List of .tre files
    :param List ali_files: List of .ali files
    :param List fasta_files: List of .fasta files

    :return file_list: List of lists containing .tre and corresponding .fasta and .ali files
    """
    file_list = []
    missmatch_list = []

    for tre_file in tre_files:
        temp = []
        length = 1
        name = re.sub(r"[^a-zA-Z0-9\.]+","_",os.path.basename(tre_file).split(".")[0])
        temp.append(tre_file)

        if ali_files:
            length += 1
            for ali_file in ali_files:
                ali = re.sub(r"[^a-zA-Z0-9\.]+","_",os.path.basename(ali_file))
                if name in ali:
                    temp.append(ali_file)

        if fasta_files:    
            length += 1  
            for fasta_file in fasta_files:
                fasta = re.sub(r"[^a-zA-Z0-9\.]+","_",os.path.basename(fasta_file))
                if name in fasta:
                    temp.append(fasta_file)
        
        if len(temp) < length:
            missmatch_list.append(temp)
        else:
            file_list.append(temp)

    return file_list, missmatch_list

def write_tsv_output(path: str,header: list, data):
    """
    Writes a tsv file from a dictionary or list.

    :param String path: Path to file that is supposed to be created
    :param List header: List of coloumn titles
    :param Dict/List data: List or dictionary of the data that should be put into .tsv file
    """

    write_header = "\t".join(header)
    with open(path, "w") as out:
        
        out.write(write_header)
        out.write("\n")

        if type(data) == list:
            for idx,i in enumerate(data):
                i = [str(x) for x in i]
                out.write("\t".join(i))
                out.write("\n")
        
        if type(data) == dict:
            for liste in [data.items()]:
                liste = [str(x) for x in liste]
                out.write("\t".join(liste))
                out.write("\n")

def write_log_file(content: list, directory: str):
    path = os.path.join(directory, "log.txt")
    with open(path, "w") as out:
        for idx,i in enumerate(content):
            if idx + 1 == len(content):
                out.write(i)
            else:
                out.write(i + "\n")
