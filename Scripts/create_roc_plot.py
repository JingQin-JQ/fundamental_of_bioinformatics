#!/usr/bin/python

# This script reads and parses your previously obtained results out of a blast output file, and a benchmark output file.
# It then creates the corresponding ROC plot, and calculates the AUC.

import argparse
import numpy
import matplotlib
matplotlib.use('AGG')
import pylab
import sys
import csv

"""
/home/student/skeleton_scripts/create_roc_plot.py -iblast /home/student/results/output.txt -ibench /home/student/results/output_classify -o /home/student/results/roc
"""

def parse_blast_results(filename):
    """
    Parse every protein pair's e-value out of a BLAST results file.
    :param filename: input file with BLAST results.
    :return: dictionary with a tuple of two UniProt IDs (key), and the corresponding e-value from BLAST (value).
    """

    blast_results = {}

    with open(filename,'r') as f:
        for line in f:
            line = line.rstrip()
            arr = line.split("\t")
    
            if len(arr) != 4:
                print("Warning: the following line does not have three elements separated by a tab:\n", line)
            elif arr[0] == arr[1]:
                print("Warning: Comparing protein to itself:", arr[0])
                
            key = (arr[0], arr[1])
            if arr[2] == "NA":    # Substitute protein pairs whose e-value is
                value = 1e6       # not available with an e-value of 1 million
            else:
                value = float(arr[2])
            blast_results[key] = value

    return blast_results


def parse_benchmark_results(filename):
    """
    Parse every protein pair's classification out of the benchmark file.
    :param filename: input file with benchmark classifacations.
    :return: dictionary with a tuple of two UniProt IDs (key), and the corresponding call (value).
    """

    benchmark_results = {}

    with open(filename,'r') as f:
        for line in f:
            line = line.rstrip()
            arr = line.split("\t")
    
            if len(arr) < 3:
                print("Warning: the following line does not have three elements separated by a tab:\n", line)
            elif arr[0] == arr[1]:
                print("Warning: Comparing protein to itself:", arr[0])

            # Benchmark classifications are symmetric, so add both possible keys:                
            key1 = (arr[0],arr[1])
            key2 = (arr[1],arr[0])
            value = arr[2]
            benchmark_results[key1] = value
            benchmark_results[key2] = value

    return benchmark_results


def integrate(x, y):
    """
    Calculate the Area Under the Curve (AUC) for a given list of coordinates
    :param x: a list of x-coordinates
    :param y: a list of y-coordinates
    :return: a float with the surface area under the curve described by x and y
    """
    auc = 0.
    last_x = x[0]
    last_y = y[0]
      
    for cur_x, cur_y in zip(x, y)[1:]:
        #########################
        ### START CODING HERE ###
        #########################
        add = (cur_x - last_x) * cur_y
  
        auc += add
        #########################
        ###  END CODING HERE  ###
        #########################
        last_x = cur_x
        last_y = cur_y
    return auc


def roc_plot(blast_evalues, benchmark_dict, png_filename, return_auc):
    """
    Draw the ROC plot for a given set of e-values and corresponding benchmark classifications.

    :param blast_evalues: the dictionary produced by parse_blast_results()
    :param benchmark_dict: the dictionary produced by parse_benchmark_results()
    """

    ### Create the lists of coordinates

    x = [0] # array of the ROC plot's x-coordinates: False Positive Rate = FP/(FP+TN)
    y = [0] # array of the ROC plot's y-coordinates: True  Positive Rate = TP/(TP+FN)

    FP = float(1e-15)
    TP = float(1e-15)
        
    # totals
    totalp = len([i for i in benchmark_dict.values() if i == 'Similar'])
    totaln = len([i for i in benchmark_dict.values() if i == 'Different'])

    last_evalue = -1
    evalues = [(v, k) for k, v in blast_evalues.items()] # List of tuples consisting of (evalue, protein_pair)
    sorted_evalues = sorted(evalues)
    for evalue, protein_pair in sorted_evalues:

        #########################
        ### START CODING HERE ###
        #########################
        # Iterate through the protein pairs, in order of ascending e-value
        # Treat every unique e-value as a homology threshold:
        #   Append one coordinate to x and y for each e-value, tracking how
        #   many "different" and "similar" pairs you've come across so far,
        #   (i.e. do not add 2 coordinates if 2 protein pairs have the same e-value.)
        # Ignore entries in the benchmark_dict classified as "ambiguous",
        #   as well as protein pairs for which you have no benchmark classification
        result_gold = benchmark_dict.get(protein_pair)
        result = ""        

        if evalue < 1:
            result = "Similar"
        else:
            result = "Different"
        if result == "Similar" and result_gold == "Similar":
            TP += float(1)
        if result == "Similar" and result_gold == "Different":
            FP += float(1)
        
                  
        x.append(FP/(totaln))
        y.append(TP/(totalp))

        #########################
        ###  END CODING HERE  ###
        #########################
        last_evalue = evalue

    # Remember ROC plots start in (0,0) and end in (1,1)!
    x = numpy.array(x) / float(x[-1]) # At this point, x[-1] = sum(benchmark_dict.values() == "different")
    y = numpy.array(y) / float(y[-1]) #                y[-1] = sum(benchmark_dict.values() == "similar")

    ### Figure out the AUC
    auc = integrate(x, y)
    
    ### Draw the plot and write it to a file
    pylab.plot(x, y)
    pylab.plot([0,1],[0,1],'--k')
    pylab.xlabel('False Positive Rate')
    pylab.ylabel('True Positive Rate')
    pylab.title('AUC = %.3f' % auc)
    pylab.savefig(png_filename)

    if return_auc:
        print auc
    return str([list(x), list(y)])

def main(blast_results_file, benchmark_results_file, png_file, return_auc):
    # Parse the input files and retrieve every protein pair's e-value and benchmark classification.
    blast_evalues = parse_blast_results(blast_results_file)

    
    benchmark_results = parse_benchmark_results(benchmark_results_file)
    
    # Draw and save the ROC plot
    print roc_plot(blast_evalues, benchmark_results, png_file, return_auc)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Draw and save a ROC plot to a file")
    parser.add_argument("-iblast","--input_blast_results", help="tab-separated BLAST results file", required=True)
    parser.add_argument("-ibench","--input_benchmark_results", help="tab-separated benchmark classification file", required=True)
    parser.add_argument("-o", "--output_png", help="output png file", required=True)
    parser.add_argument("-auc", "--return_auc", dest="return_auc", action="store_true", help="If flagged, return AUC")

    args = parser.parse_args()
    
    blast_file = args.input_blast_results
    benchmark_file = args.input_benchmark_results
    png_file = args.output_png
    return_auc = args.return_auc 

    main(blast_file,benchmark_file, png_file, return_auc)
