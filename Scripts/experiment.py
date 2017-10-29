#!/usr/bin/python

"""
experiment.py

Runs the blast vs psiblast experiment:
- Runs run_local_blast.py multiple times:
    - For blast once
    - For psiblast multiple times with varying parameters
- Checks the results against the gold standard:
    - Produces ROC plots for the blast and psiblast results
"""

# import packages
import os
import sys
import subprocess
import numpy
import time
import csv
import matplotlib
matplotlib.use('AGG')
import pylab

# get arguments for psiblast parameters: number of iterations and evalue inclusion threshold
if len(sys.argv) == 4 and len(sys.argv[1].split(',')) == 3 and len(sys.argv[2].split(',')) == 3:
    num_iter_start = int(sys.argv[1].split(',')[0])
    num_iter_end = int(sys.argv[1].split(',')[1])
    num_iter_step = int(sys.argv[1].split(',')[2])
    e_start = float(sys.argv[2].split(',')[0])
    e_end = float(sys.argv[2].split(',')[1])
    e_step = float(sys.argv[2].split(',')[2])
    evalue_incl_thresh = numpy.arange(e_start,e_end + e_step,e_step)
    gold_method = sys.argv[3]    
    
else:
    sys.exit("Please run with three arguments: [start,end] [start,end,step] [fold/all]\nWhere [start,end] are the start and end integers for the number of iterations and [start,end,step] are the start, end and step values for the evalue inclusion threshold range. [fold/all] is the golden standard method used.")

## BLAST
blast_command = "/home/student/skeleton_scripts/run_local_blast.py \
       -ids /home/student/data/uniprot_id_list.txt \
       -q /home/student/queries \
       -db /home/student/db/complete.fasta \
       -o /home/student/results/blast/%s/output_blast.txt \
       -opng /home/student/results/blast/%s/DistributionEValue_blast.png" % (gold_method,gold_method)

# run blast benchmark
print "Running blast benchmark..."
subprocess.call(blast_command, shell = True)
blast_output = subprocess.check_output("/home/student/skeleton_scripts/create_roc_plot.py -iblast /home/student/results/blast/%s/output_blast.txt -ibench /home/student/results/output_classify_%s -o /home/student/results/blast/%s/roc_blast -auc" % (gold_method,gold_method,gold_method), shell = True)

# save xy_blast roc points to variable
xy_blast = eval(blast_output.split('\n')[1])
for i in range(len(xy_blast)):
    xy_blast[i] = numpy.array(xy_blast[i])
    

print "BLAST AUC: %s\n" % str(blast_output.split('\n')[0]) + "Finished blast benchmark...\n"

## PSIBLAST
# counter for status print
count = 0

# result_auc list to store auc results
result_auc = []

# list to store all xy coordinates of returned roc values
xy_psiblast = []

# run psiblast benchmark with varying parameters
print "Running psiblast benchmark with parameter optimalization..."
for n in range(num_iter_start,num_iter_end + 1, num_iter_step):
    for i in evalue_incl_thresh:

        # counter for status print
        count += 1
        
        # print status
        sys.stdout.write("\r" + "Running psiblast %s of %s..." % (count, len(evalue_incl_thresh) * len(range(num_iter_start,num_iter_end + 1, num_iter_step))))
        sys.stdout.flush()
           
        # run psiblast with the parameters n and i  
        psiblast_command = "/home/student/skeleton_scripts/run_local_blast.py \
           -ids /home/student/data/uniprot_id_list.txt \
           -q /home/student/queries \
           -db /home/student/db/complete.fasta \
           -o /home/student/results/psiblast/%s/output_psiblast_evalue_incl_thresh_%s_numiter_%s.txt \
           -opng /home/student/results/psiblast/%s/DistributionEValue_psiblast_evalue_incl_thresh_%s_numiter_%s.png \
           -psi \
           -n %s \
           -incl %s" % (gold_method,str(i).replace('.',','),n,gold_method,str(i).replace('.',','),n,n,i)
        subprocess.call(psiblast_command, shell = True)
        
        # run the create_roc_plot script and save auc result
        psiblast_output = subprocess.check_output("/home/student/skeleton_scripts/create_roc_plot.py -iblast /home/student/results/psiblast/%s/output_psiblast_evalue_incl_thresh_%s_numiter_%s.txt -ibench /home/student/results/output_classify_%s -o /home/student/results/psiblast/%s/roc_psiblast_evalue_incl_thresh_%s_numiter_%s -auc" % (gold_method,str(i).replace('.',','),n,gold_method,gold_method,str(i).replace('.',','),n) ,shell = True)
        
        # add auc result to result_auc dict and save xy_psiblast roc points
        result_auc.append([i,n,float(psiblast_output.split('\n')[0])])
        xy_psiblast.append(eval(psiblast_output.split('\n')[1]))

print "\nFinished psiblast benchmark..."

# get maximum auc value from result_auc list
maximum = max([sublist[-1] for sublist in result_auc])

# find parameters in result_auc for max auc
maximum_index = 0
for i in result_auc:
    if i[2] == maximum:
        param = i
        break
    maximum_index += 1

# save result_auc to .csv
with open("/home/student/results/psiblast/%s/psiblast_result_%s_%s.csv" % (gold_method,sys.argv[1],sys.argv[2]), "wb") as f:
    writer = csv.writer(f)
    writer.writerow(["Blast AUC: %s" % str(blast_output.split('\n')[0]).rstrip()])
    writer.writerow(["Max Psiblast AUC: %s at e value inclusion shreshold = %s and number of iterations = %s" % (str(maximum), str(param[0]), str(param[1]))])
    writer.writerow(["\nPsiblast parameter iteration:\n"])
    writer.writerow(["e value inclusion shreshold","Number of iterations","AUC"])
    writer.writerows(result_auc)

# make plot for blast
pylab.plot(xy_blast[0],xy_blast[1], label = 'BLAST, AUC = %.3f' % float(blast_output.split('\n')[0]),linewidth=2)

# convert xy_psiblast to numpy.array with floats
for i in range(len(xy_psiblast)):
    for xy in range(len(xy_psiblast[i])): 
        xy_psiblast[i][xy] = numpy.array(xy_psiblast[i][xy])

# make plot for psiblast
pylab.plot(xy_psiblast[maximum_index][0],xy_psiblast[maximum_index][1], label = 'PSI-BLAST, AUC = %.3f, num_iter = %s, e_incl_thresh = %s' % (maximum,str(result_auc[maximum_index][1]),str(result_auc[maximum_index][0])),linewidth=2)

# plot diagonal
pylab.plot([0,1],[0,1],'--k')

# add labels, title and legend
pylab.xlabel('False Positive Rate')
pylab.ylabel('True Positive Rate')
pylab.title('BLAST vs PSI-BLAST\nCompared SCOP data term(s): %s' % gold_method)
pylab.legend(loc='lower right', prop={'size': 7})
# save plot
pylab.savefig('/home/student/results/roc_combined_%s.png' % gold_method)

# print result
print "\n-----------------------------------------------------"
print "----------------------PSIBLAST-----------------------"
print "-----------Parameter optimalization result-----------\n"
print "e value inclusion shreshold: " + str(param[0]) + "\n" + "Number of iterations: " + str(param[1]) + "\nResulting AUC: " + str(maximum)
print "Results saved in ~/results/blast/%s and ~/results/psiblast/%s!" % (gold_method,gold_method)
