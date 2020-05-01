"""========================================================
	  		FOF FASTGWA-MLM PIPELINE V1.0
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
import cgatcore.experiment as E
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


@originate(PARAMS["inp_title"]) 
@mkdir(PARAMS["GRM_dir"])         
def generate_GRM(title):	


	'''
	Function to generate Genetic Relation Matrix from the
	FASTGWA-MLM analysis.
	'''	


	GRM_dir = PARAMS["GRM_dir"]


	statement = '''

	~/Desktop/gcta/gcta64 --bfile %(title)s
	--make-grm --sparse-cutoff 0.05 --threads 10
	--out %(GRM_dir)s/%(title)s

				'''

	P.run(statement)



##################################################
# Run mlm-fastgwa
##################################################



@originate(PARAMS["inp_title"]) 
@follows(generate_GRM)
@mkdir(PARAMS["asso_dir"])
def mlm_gwa(title):

	'''
	Function to run fastGWA on the plink files and GRM.
	'''	


	GRM_dir = PARAMS["GRM_dir"]
	asso_dir = PARAMS["asso_dir"]
	pheno = PARAMS["pheno"]


	statement = ''' 

	~/Desktop/gcta/gcta64 --bfile %(title)s
	--grm-sparse %(GRM_dir)s/%(title)s --fastGWA-mlm 
	--pheno %(pheno)s --threads 10 
	--out %(asso_dir)s/%(title)s

				'''

	P.run(statement)



#Path of the output of fastGWA
gwas_file = os.path.abspath(
        os.path.join(PARAMS["asso_dir"] + "/*.fastGWA"))



@follows(mlm_gwa)
@transform([gwas_file], suffix('.fastGWA'), '.tsv')
def reformat(infile, outfile):

	'''
	Function to reformat the fastGWA output so it is readable
	into Hail for Manhattan and QQ plotting. Currently missing
	code to reformat locus (Manhattan plot not currently 
	working).
	'''	

	df = pd.read_csv(infile, sep='\t')
	
	df_list = df['SNP'].tolist()
	df_list = [w.replace('rs', '') for w in df_list]
	df['SNP'] = pd.Series(df_list)

	df["locus"] = df["CHR"].astype(str) + ':' + df["SNP"].astype(str)

	p_list = df['P'].tolist()
	bp_list = []

	for i in p_list:
		t = np.format_float_scientific(np.float32(i))
		bp_list.append(t)

  
	df['P'] = pd.Series(bp_list)

	df.to_csv(outfile, index=None, sep ='\t')


	pass



#Path of the newly formatted fastGWA output
gwas_f_file = os.path.abspath(
        os.path.join(PARAMS["asso_dir"] + "/*.tsv"))


@follows(reformat)
@transform([gwas_f_file], suffix('.tsv'), '.png')
def QQ_plot(infile, outfile):

	'''
	Function to perform QQ_plot of p-vals of fastGWA.
	'''	

	mlm = hl.import_table(infile,
                              types={'P': hl.tfloat64})

	p = hl.plot.qq(mlm.P, title="QQ_plot")

	export_png(p, filename=outfile)


	pass



@follows(QQ_plot)
def full():
	'''
	Dummy function.
	'''	
    pass



#Pass arguments
def main(argv=None):
    if argv is None:
        argv = sys.argv
    P.main(argv)


if __name__ == "__main__":
    sys.exit(P.main(sys.argv))






