"""========================================================
	  		FOF FASTGWA-MLM PIPELINE V1.1
===========================================================

Overview
========
The pipeline performs the following:
   * GRM generation from plink files
   * Mixed linear model fastGWA analysis
   * Manhattan and QQ plot

Author: Felix O'Farell
Date: April

   
   
Usage
=====



Configuration
-------------
The pipeline requires a pipeline.yml configuration file. This is located
within the fast-mlm directory.
Input
-----
The title of the plink files. E.g. test.bim, test.fam, test.bam --> title = "test"

Output
-----
* fastGWA tsv file
* QQ plot
* Manhattan plot


Code
====
"""

import cgatcore.experiment as E
from ruffus import *
from ruffus.combinatorics import *
from cgatcore import pipeline as P
import cgatcore.iotools as iotools
import os
import sys
import pandas as pd
import numpy as np
import hail as hl
from bokeh.io import export_png




#get parameters from YML file
PARAMS = P.get_parameters(
    ["%s/pipeline.yml" % os.path.splitext(__file__)[0],
    "../pipeline.yml",
    "pipeline.yml"])


##################################################
# Generate GRM from raw plink files
##################################################


c_dir = PARAMS["inp_file"]
c_num = PARAMS["chr_number"]
c_title = PARAMS["inp_title"]

chrx = []

#build list of inputs for GRM and fastGWA functions
for i in range(1,c_num):
	name = c_dir + '/' + c_title +  str(i) + '/' +  c_title + str(i)
	chrx.append(name)


##################################################
# Generate GRM from raw plink files
##################################################


@originate(chrx) 
@mkdir(PARAMS["GRM_dir"])          
def generate_GRM(title):	

	'''
	Function to generate Genetic Relation Matrix from the
	FASTGWA-MLM analysis.
	'''	

	GRM_dir = PARAMS["GRM_dir"]

	num = title[(len(title)-2):len(title)]

	num = num.replace('r','')


	statement = '''

	~/Desktop/gcta/gcta64 --bfile %(title)s
	--make-grm --sparse-cutoff 0.05 --threads 10
	--out %(GRM_dir)s/chr%(num)s

				'''


	P.run(statement)


##################################################
# Run mlm-fastgwa
##################################################



@follows(generate_GRM)
@originate(chrx) 
@mkdir(PARAMS["asso_dir"])
def mlm_gwa(title):

	'''
	Function to run fastGWA on the plink files and GRM.
	'''	

	GRM_dir = PARAMS["GRM_dir"]

	asso_dir = PARAMS["asso_dir"]

	pheno = PARAMS["pheno"]

	num = title[(len(title)-2):len(title)]

	num = num.replace('r','')


	statement = ''' 

	~/Desktop/gcta/gcta64 --bfile %(title)s
	--grm-sparse %(GRM_dir)s/chr%(num)s --fastGWA-mlm 
	--pheno %(title)s.fam --threads 10 
	--out %(asso_dir)s/chr%(num)s

				'''


	P.run(statement)



##################################################
# Create master fastGWA
##################################################


@follows(mlm_gwa)
@mkdir('masterGWA')
def master_maker():
	

	'''
	Function to copy Chr.fastGWA as masterGWA file 
	for concat function.
	'''	


	statement = '''


	cp asso/chr1.fastGWA masterGWA/master.fastGWA


				'''
	

	P.run(statement)


##################################################
# Concat the individual fastGWA's
##################################################


gwas_files = os.path.abspath(
       os.path.join(PARAMS["asso_dir"] + "/*.fastGWA"))

#Need to stipulate 1 thread for this function
@follows(master_maker) 
@transform([gwas_files], suffix('.fastGWA'), '.tsv')
def concat(infile, outfile):

	'''
	Function to append the previous .fastGWA summary
	file to the next. This outputs a masterfastGWA 
	summary file which contains SNPs across entire
	genome.
	'''	


	statement = '''


	Rscript concat.r  masterGWA/master.fastGWA  %(infile)s 


				'''
	P.run(statement)


##################################################
# Reformat the master summary file for Hail
##################################################


#masterGWA file path
masterGWA_f = os.path.abspath(
       os.path.join('masterGWA' + "/*.fastGWA"))


@follows(concat)
@transform([masterGWA_f], suffix('.fastGWA'), '_h.tsv')
def hail_format(infile, outfile):


	'''
	Function to reformat the masterfastGWA output so it is readable
	into Hail for Manhattan and QQ plotting. Currently missing
	code to reformat locus (Manhattan plot not currently 
	working).
	'''	

	statement = '''


    Rscript hail_prep.r %(infile)s %(outfile)s


                '''


	P.run(statement)


#Path of the newly formatted fastGWA output


##################################################
# Create QQ plot  
##################################################

gwas_m_file = os.path.abspath(
        os.path.join('masterGWA' + "/*.fastGWA"))

@follows(hail_format)
@transform([gwas_m_file], suffix('.fastGWA'), '_QQ.png')
def QQ_plot(infile, outfile):


	'''
	Function to perform QQ_plot of p-vals of masterfastGWA.
	'''	

	mlm = hl.import_table(infile,
                        types={'P': hl.tfloat64,
                      'locus': hl.tlocus(reference_genome = 'GRCh37')})


	p = hl.plot.qq(mlm.P, title='QQ_plot')

	export_png(p, filename=outfile)

	

	pass

##################################################
# Create Manhattan plot
##################################################


@follows(QQ_plot)
@transform([gwas_m_file], suffix('.fastGWA'), '_Man.png')
def man_plot(infile, outfile):


	'''
	Function to perform man_plot of p-vals of masterfastGWA.
	'''	


	mlm = hl.import_table(infile,
                        types={'P': hl.tfloat64,
                      'locus': hl.tlocus(reference_genome = 'GRCh37')})


	f = hl.plot.manhattan(mlm.P, title='Man_Plot')

	export_png(f, filename=outfile)


	pass


##################################################
# Clean up function 
##################################################


@follows(man_plot)
@mkdir('logs')
def clean_up():


	'''
	Function to move log files into a log directory 
	'''

	statement = '''


	mv *.sh *.times logs


				'''
	P.run(statement)


##################################################
# Dummy function 
##################################################


@follows(clean_up)
def full():
	

    pass


##################################################
# Pass aguments 
##################################################


#Pass arguments
@follows(full)
def main(argv=None):
    if argv is None:
        argv = sys.argv
    P.main(argv)


if __name__ == "__main__":
    sys.exit(P.main(sys.argv))
