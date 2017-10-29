#!/usr/bin/python

import itertools
import argparse
import re
import sys

"""./classify_scop_skeleton.py -o /home/student/results/output_classify -ids /home/student/data/uniprot_id_list.txt -pdb /home/student/data/PDB_ID_lookup.tab -s /home/student/data/dir.cla.scop.txt
"""

def retrieve_scop_data(scop_file):
    """
    Reads a databaset file from SCOP and returns a dictionary with the protein IDS mapping.
    :param scop_file: database file containing mapping of PDB's to SCOP ID's.
    :return: dictionary containing the mapping of PDBs (keys) and their corresponding SCOP data (values).
    """

    scop_data = {}

    ##########################
    ### START CODING HERE ####
    ##########################
    # You can parse SCOP data in various ways. E.g. you can use dictionary of dictionaries
    # {proteinID: {"class": class, "fold": fold, "superfamily": superfamily, 'family': family}}
    scop_file = open('/home/student/data/dir.cla.scop.txt')

    # Skip the first 4 lines
    for i in xrange(4):
        scop_file.next()

    # iterate all data lines in the scop data file
    for line in scop_file:
        protein_data = line.strip()

        # serparate the elements in the line
        proteinID = protein_data.split("\t")[1]
        info = protein_data.split("\t")[5]
        info_list= info.split(",")

        # extract profile numbers
        for i in range(len(info_list)):
            info_list[i] = info_list[i][3:]

        # add each element to the dictionary
        l= ["class","fold","superfamily","family"]
        scop_data[proteinID] = {} 
        scop_data[proteinID]=dict(zip(l,info_list))

    # close the scop file
    scop_file.close()
    ########################
    ### END CODING HERE ####
    ########################

    return scop_data


def compute_similarity_score(prot1_scop, prot2_scop):
    """
    Computes the score for two proteins on the basis of the data from SCOP database.
    :param prot1_scop: data for protein 1 from SCOP database.
    :param prot2_scop: data for protein 2 from SCOP database.
    :return: similarity score (float)
    """
    ##########################
    ### START CODING HERE ####
    ##########################
    """
    Calculate similarity score (css):
    Compare each pdb id data form protein 1 to each pdb id data from protein 2
    If all data elements are the same for a comparrison (family, superfamily, fold etc), add 1 to the score
    Divide score by the number of domains in protein 1 * the number of domains in protein 2
    """
    score1 = float(0)
    score2 = float(0)
    score_part = 0
    score = float(0)

    # if no data for a protein: return 0
    if len(prot1_scop) == 0 or len(prot2_scop) == 0:
        return score
    
    # compare pdb's for each protein and add 1 to score if they match, do this forward and backward
    
    for i in prot1_scop:
        for j in prot2_scop:
            if i['fold'] == j['fold']:
                score_part = 1
        if score_part == 1:
            score1 += 1
        score_part = 0
    score1 /= len(prot1_scop)
    
    for i in prot2_scop:
        for j in prot1_scop:
            if i['fold'] == j['fold']:
                score_part = 1
        if score_part == 1:
            score2 += 1   
        score_part = 0 
    score2 /= len(prot2_scop) 
        
    # sum the partial scores and find the average     
    score = score1 + score2
    score /= 2

    return score
    ########################
    ### END CODING HERE ####
    ########################


def check_similarity_for_protein_pair(prot1_scop, prot2_scop):
    """
    Returns the similarity score between two proteins.
    :param prot1_scop: data for protein 1 from SCOP database.
    :param prot2_scop: data for protein 2 from SCOP database.
    :param pair: a tuple with the UniProt IDs of the two proteins to compare.
    :return: "different", "similar" or "ambiguous".
    """
    ##########################
    ### START CODING HERE ####
    ##########################
    pass
    ########################
    ### END CODING HERE ####
    ########################

# If you will use the numeric score for SCOP (similar to GO), you may want to use check_similarity_for_protein_pair
# with other arguments. See the example below.
# def check_similarity_for_protein_pair(score, threshold):
#    pass

def generate_all_possible_protein_pairs(protein_ids):
    """
    Returns a list containing all unique protein pairs.
    :param protein_ids: list of all proteins IDs.
    :return: list of possible unique protein pairs.
    """
    pairs = list()
    ##########################
    ### START CODING HERE ####
    ##########################
    for i in itertools.combinations(protein_ids,2):
        pairs.append(i)

    ########################
    ### END CODING HERE ####
    ########################
    return pairs


def assign_homology(scop_dict, protein_ids_pdbs, pairs):
    """
    Computes the similarity score between all protein pairs from the list, and decides if two proteins are homologs
    (different, ambiguous or similar).
    :param scop_dict: dictionary containing the mapping of PDBs (keys) and their corresponding SCOP data (values).
    :param protein_ids_pdbs: dictionary with UniprotID as key and PDB ID as a value.
    :param pairs: list of all possible unique protein pairs.
    :return: dictionary with UniProt ID (key), similarity(different, ambiguous or similar).
    """
    scop_homology = {}

    ##########################
    ### START CODING HERE ####
    ##########################
    # You should remember to take care about the proteins that are not in the SCOP database.

    for j in pairs:

        # retrieve pdb's for each protein
        prot1_pdb = protein_ids_pdbs.get(j[0])
        prot2_pdb = protein_ids_pdbs.get(j[1])

        prot1_scop = []        
        prot2_scop = []

        # add all pdb id data entries to the data list
        for i in prot1_pdb:
            prot1_scop.append(scop_dict.get(i.lower()))
        for i in prot2_pdb:
            prot2_scop.append(scop_dict.get(i.lower()))
        
        # flag proteins without SCOP data (serperate loops for each protein to enable using the "continue" statement
        any_data = False
        for domain in prot1_scop:
            if domain != None:
                any_data = True
        if any_data == False:
            scop_homology[(j[0],j[1])] = "Ambiguous"
            continue

        any_data = False
        for domain in prot2_scop:
            if domain != None:
                any_data = True
        if any_data == False:
            scop_homology[(j[0],j[1])] = "Ambiguous"
            continue

        # take out None entries (pdb's without data)
        prot1_scop = filter(None, prot1_scop)
        prot2_scop = filter(None, prot2_scop)

        # get the similarity score for the pdb's
        score = compute_similarity_score(prot1_scop, prot2_scop)

        # add the score for the proteins to the dictionary
        scop_homology[(j[0],j[1])] = score
           
    ########################
    ### END CODING HERE ####
    ########################

    return scop_homology


def write_results(filename, scop_homology):
    """
    Writes in an output file the all of the protein pairs and their similarity/dissimilarity.
    :param output_file: the name of the output file.
    :param scop_homology: dictionary (keys: protein pairs as tuples; values: one of the value - different/similar/ambiguous)
    """
    with open(filename, "w") as f:
        for (p1, p2), value in scop_homology.iteritems():
            if value == "Ambiguous":
                f.write("\t".join([str(p1), str(p2), "Ambiguous"]) + "\n")
            elif float(value) > 0.5:
                f.write("\t".join([str(p1), str(p2), "Similar"]) + "\n")
            else:
                f.write("\t".join([str(p1), str(p2), "Different"]) + "\n")


def read_protein_ids_file(filename):
    """
    Returns the list of UniProt IDs contained in the input file.
    :param filename:
    :return: list of UniProt IDs in the input file.
    """
    ##########################
    ### START CODING HERE ####
    ##########################
    protein_ids = []
    uniprot_ids_file = open('/home/student/data/uniprot_id_list.txt')
    for line in uniprot_ids_file:
        query = line.strip()
        protein_ids.append(query)

    uniprot_ids_file.close()
    #######################
    ### END CODING HERE ###
    #######################
    return protein_ids


def read_lookup_table(filename):
    """
    Reads the specified file and returns the dictionary with UniprotID as key and PDB ID as a value.
    :param filename: file with the mapping between Uniprot ids and PDB ids.
    :return: dictionary with UniprotID as key and PDB ID as a value.
    """
    ##########################
    ### START CODING HERE ####
    ##########################
    pdb_id_file = open('/home/student/data/PDB_ID_lookup.tab')
    protid_pdb = {}

    # add each uniprot and pdb to a dict
    for line in pdb_id_file:
        line_strip = line.strip()
        protid_key = line_strip.split("\t")[0]
        protid_value = line_strip.split("\t")[1]

        if protid_pdb.has_key(protid_key):
            protid_pdb[protid_key].append(protid_value)
        else:            
            protid_pdb[protid_key]= [protid_value]

    pdb_id_file.close()
    #######################
    ### END CODING HERE ###
    #######################
    return protid_pdb


def main(input_file, output_file, pdb_id_file, scop_file):
    ##########################
    ### START CODING HERE ####
    ##########################
    
    # get uniprot id's
    protein_ids = read_protein_ids_file(input_file)
    
    # get uniprot to pdb conversion table
    protid_pdb = read_lookup_table(pdb_id_file)

    # generate pairs
    pairs = generate_all_possible_protein_pairs(protein_ids)

    # get scop data
    scop_data = retrieve_scop_data(scop_file)

    # get homologies
    scop_homology = assign_homology(scop_data, protid_pdb, pairs)

    # write resuts
    write_results(output_file, scop_homology)
    
    #######################
    ### END CODING HERE ###
    #######################


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='The script retrieves data from SCOP database (from local file)'
                                                 ' and provides an output file'
                                                 ' with the strings "ProteinID   ProteinID   similarity", where'
                                                 ' similarity is a string with one of the values from '
                                                 ' different/ambiguous/similar.')

    parser.add_argument("-o", "--output_file", help="Output file name")
    parser.add_argument("-ids", "--protein_ids_file", help="File with the protein Uniprot IDs")
    parser.add_argument("-pdb", "--pdb_id_file", help="File with the mapping between Uniprot ids and PDB ids")
    parser.add_argument("-s", "--scop_file", help="SCOP database file")

    args = parser.parse_args()

    input_file = args.protein_ids_file
    output_file = args.output_file
    pdb_id_file = args.pdb_id_file
    scop_file = args.scop_file

    main(input_file, output_file, pdb_id_file, scop_file)
