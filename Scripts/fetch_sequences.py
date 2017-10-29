#!/usr/bin/python
# -*- coding: utf-8 -*-

# This script downloads all the sequences in the list provided by the user, and puts them into a single file to be used
# as a database.
# It also stores the fetched fasta sequences into individual files which will be used as queries for (PSI-)BLAST search.

import urllib2
import argparse
import os.path


def fetch_one_fasta(uniprot_id):
    """
    Fetch the fasta formatted sequence of the uniProtID supplied in the argument.
    :param uniProtID: id of the protein to fetch in uniprot database.
    :return: sequence data in fasta format.
    """
    url = 'http://www.uniprot.org/uniprot/' + uniprot_id + '.fasta'

    fh = urllib2.urlopen(url)
    result = fh.read()
    fh.close()
    return result


def check_query_folder(query_folder):
    if not os.path.exists(query_folder):
        return False
    return True


def fetch_all_sequences(query_folder, uniprot_ids, database_file):
    """
    Fetch the fasta formatted sequence for each uniProtID in the list.
    :param query_folder: folder containing fasta files.
    :param uniprot_ids: File with the list of protein IDs
    :param database_file: database files to save the results of the fetching.
    """
    if not query_folder.endswith("/"):
        query_folder += "/"

    print("Processing the list of ids...")
	
	# # for the .fasta file containing all fasta sequences
    f2 = open("/home/student/db/complete.fasta",'a')
	
    for line in uniprot_ids:
		
        fasta = fetch_one_fasta(line.rstrip())
        
		# for all individual .fasta files; add to query folder
        f = open('/home/student/queries/%s' % line.rstrip() + '.fasta','w')
        f.write(fasta)
        f.close()		

		# for the .fasta file containing all fasta sequences
        f2.write(fasta)

	# close complete fasta file
    f2.close()

    print("Processing finished.")
    uniprot_ids.close()
    database_file.close()


def main(query_folder, uniprot_ids, database_file):
    if check_query_folder(query_folder):
        fetch_all_sequences(query_folder, uniprot_ids, database_file)
    else:
        print("Query folder does not exist. "
              "Please make sure that the specified folder exists before you run this script.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='This script downloads all the sequences in the list provided'
                                                 'by the user, and puts them into a single file to be used'
                                                 'as a database. It also stores the fetched fasta sequences into'
                                                 'individual files which will be used as queries for (PSI-)BLAST search')

    parser.add_argument("-ids", "--uniprot_ids", type=argparse.FileType("r"),
                        help="File with the list of protein IDs", required=True)
    parser.add_argument("-db", "--dbfile", type=argparse.FileType("w"),
                        help="Output file containing all fetched sequences", required=True)
    parser.add_argument("-q", "--qfolder", help="Output folder for your queries", required=True)
    args = parser.parse_args()

    query_folder = args.qfolder
    uniprot_ids = args.uniprot_ids
    database_file = args.dbfile

    main(query_folder, uniprot_ids, database_file)
