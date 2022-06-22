import argparse
from importlib.resources import path
import logging
from pathlib import Path
import subprocess
from Bio import SeqIO
import numpy as np
from sympy import re
from tqdm import tqdm
import math
from collections import Counter
import pandas as pd
import os
import pysam
tqdm.pandas()


parser = argparse.ArgumentParser(description='Process files.')
parser.add_argument("-r", "--reference", type = Path)
parser.add_argument("--bathmet2", type = Path)
parser.add_argument("--bed",type = Path)
parser.add_argument("--bismark",type = Path)
parser.add_argument("--deepsignalplant",type = Path)
parser.add_argument("-d", "--debug", dest="loglevel", help="For debugging in log.", action="store_const", const=logging.DEBUG, default=logging.INFO)

args = parser.parse_args()

#decorateur
def print_info(func):
    def decorateur(*args):
            logging.info('Starting ' + func.__name__+' processing')
            func(*args)
            logging.info('Finished ' + func.__name__+' processing')
    return decorateur

def transfoContext(row):
    #https://stackoverflow.com/questions/26886653/pandas-create-new-column-based-on-values-from-other-columns-apply-a-function-o
    if row['Contexte2'][3]=='G':
        return 'CG'
    elif row['Contexte2'][4]=='G':
        return 'CHG'
    else : return 'CHH'

def GetBed(infile,outfile):
    nrow, _ = subprocess.check_output(["wc", "-l", str(infile)]).decode().split()
    logging.info("total rows: "+nrow)
    res = pd.concat([chunk for chunk in tqdm(pd.read_csv(infile, sep = "\t", usecols=[0,1,5,9,10], header=0, names=['Chr','Pos','Brin','Couverture','Pourcentage_de_methylation'], chunksize=5000), desc='Loading Bismark data', total=int(nrow)/5000)])
    res['Reads_methyles']= round(res['Pourcentage_de_methylation'] * res['Couverture']/100)
    res['Reads_non_methyles']= round(res['Couverture']-res['Reads_methyles'])
    res['Contexte']='CG'
    res['Pos_Globale']= res['Pos']+1
    res.iloc[:, [0,8,2,5,6,4,3,7]].to_csv(outfile, sep='\t', index=False)

def getBathMet2(infile,outfile):
    # bathmet :#chromsome	loci	strand	context	C_count	CT_count	methRatio	eff_CT_count	rev_G_count	rev_GA_count	MethContext	5context
    nrow, _ = subprocess.check_output(["wc", "-l", str(infile)]).decode().split()
    logging.info("total rows: "+nrow)
    res = pd.concat([chunk for chunk in tqdm(pd.read_csv(infile, sep = "\t", usecols=[0,1,2,3,4,5,6], header=0, names=['Chr','Pos_Globale','Brin','Contexte','Reads_methyles','Reads_non_methyles','Ratio_de_meth'], chunksize=5000), desc='Loading Bathmet2 data', total=int(nrow)/5000)])
    res['Couverture']= res['Reads_methyles'] + res['Reads_non_methyles'] 
    res['Pourcentage_de_methylation']= round(res ['Reads_methyles'] / res['Couverture']*100, ndigits=2)
    res = res.fillna(np.NaN) # car vide si division par 0
    return res.iloc[:, [0,1,2,4,5,8,7,3]].to_csv(outfile, sep='\t')

def getBismark(infile, outfile):
    # bismark : "chr", "pos", "strand", "methylated", "unmethylated", "context", "trinucleotide"
    nrow, _ = subprocess.check_output(["wc", "-l", str(infile)]).decode().split()
    logging.info("total rows: "+nrow)
    res = pd.concat([chunk for chunk in tqdm(pd.read_csv(infile, sep = "\t", usecols=[0,1,2,3,4,5], header=0, names=['Chr','Pos_Globale','Brin','Reads_methyles','Reads_non_methyles','Contexte'], chunksize=5000), desc='Loading Bismark data', total=int(nrow)/5000)])
    #res = pd.read_csv(infile, sep = "\t", usecols=[0,1,2,3,4,5], header=0, names=['Chr','Pos_Globale','Brin','Reads_methyles','Reads_non_methyles','Contexte'])
    res['Couverture']= res['Reads_methyles'] + res['Reads_non_methyles'] 
    res['Pourcentage_de_methylation']= round(res ['Reads_methyles'] / res['Couverture']*100, ndigits=2)
    res = res.fillna(np.NaN)
    logging.info('Bismark traite saving in '+str(outfile))
    res.iloc[:, [0,1,2,3,4,7,6,5]].to_csv(outfile, sep='\t')
    logging.info('Bismark traite saved')

def getDeepSignal(infile, outfile):
    nrow, _ = subprocess.check_output(["wc", "-l", str(infile)]).decode().split()
    logging.info("total rows: "+nrow)
    res = pd.concat([chunk for chunk in tqdm(pd.read_csv(infile, sep = "\t", usecols=[0,1,2,6,7,8,10], header=0, names=['Chr','Pos','Brin','Reads_methyles','Reads_non_methyles','Couverture','Contexte2'], chunksize=5000), desc='Loading DeepSignal data', total=int(nrow)/5000)])
    res['Pourcentage_de_methylation']= round(res ['Reads_methyles'] / res['Couverture']*100, ndigits=2)
    res['Contexte']= res.progress_apply(lambda row: transfoContext(row), axis=1)
    res['Pos_Globale']= res['Pos']+1 #car pas même base
    res = res.fillna(np.NaN)
    logging.info('DeepSignal traite saving in '+str(outfile))
    res.iloc[:, [0,9,2,3,4,7,5,8]].to_csv(outfile, sep='\t')
    logging.info('DeepSignal traite saved')


def getMerged(infile1, infile2, outfile):
    nrow, _ = subprocess.check_output(["wc", "-l", str(infile1)]).decode().split()
    logging.info("total rows of the reference: "+nrow)
    ref = pd.concat([chunk for chunk in tqdm(pd.read_csv(infile1, sep='\t', dtype={'Pos_Globale': str}, chunksize=5000), desc='Loading reference data', total=int(nrow)/5000)])
    nrow2, _ = subprocess.check_output(["wc", "-l", str(infile2)]).decode().split()
    logging.info("total rows of the sequencing file: "+nrow2)
    autre = pd.concat([chunk for chunk in tqdm(pd.read_csv(infile2, index_col=0, sep='\t', dtype={'Pos_Globale': str}, chunksize=5000), desc='Loading sequencing file data', total=int(nrow)/5000)])
    a=ref.merge(autre, how='left', on=['Chr','Contexte','Pos_Globale'])
    logging.info('Merged file saving in '+str(outfile))
    a.to_csv(outfile, sep='\t',index=False)
    logging.info('Merged file saved')

def getMergedPourcentage(infileref,infile2,outfile):  
    ref = pd.read_csv(infileref, sep='\t', usecols=['Chr','Pos_Globale','ID'])
    autre= pd.read_csv(infile2, sep='\t', usecols=['Chr','Pos_Globale','ID','Pourcentage_de_methylation'])
    a=ref.merge(autre, how='left', on=['Chr','Pos_Globale','ID'])
    a.to_csv(outfile, sep='\t', index=False)

@print_info
def analyseBathMet2(input, ref):
    getBathMet2(input,input.parent(input.stem+"Traite.txt"))
    getMerged(ref, input.parent / (input.stem+"Traite.txt"), input.parent(input.stem+"Merged.txt"))

@print_info
def analyseBismark(input, ref):
    getBismark(input,input.parent / (input.stem + "Traite.txt"))
    getMerged(ref, input.parent /(input.stem+"Traite.txt"), input.parent /(input.stem+"Merged.txt"))

@print_info
def analyseBed(input, ref):
    GetBed(input,input.parent /(input.stem+"Traite.txt"))
    getMerged(ref, input.parent /(input.stem+"Traite.txt"), input.parent /(input.stem+"Merged.txt"))

@print_info
def analyseDeepSignalPlant(input, ref):
    getDeepSignal(input, input.parent /(input.stem+"Traite.txt"))
    getMerged(ref, input.parent /(input.stem+"Traite.txt"), input.parent /(input.stem+"Merged.txt"))


def check_args(args):
    if args.bathmet2 and not args.bathmet2.exists():
        logging.warning("BathMet2 input does not exist. Please give a valid file.")
        exit()   
    if args.bed and not args.bed.exists():
        logging.warning("Bed input does not exist. Please give a valid file.")
        exit()  
    if args.bismark and not args.bismark.exists():
        logging.warning("Bismark input does not exist. Please give a valid file.")  
        exit()
    if args.deepsignalplant and not args.deepsignalplant.exists():
        logging.warning("DeepSignalPlant input does not exist. Please give a valid file.")
        exit() 



if __name__ == '__main__':
    check_args(args)
    logging.basicConfig(filename='myappdata.log', filemode='w', format='%(asctime)s %(message)s', datefmt='%d/%m/%Y %H:%M:%S', level=args.loglevel)
    logging.info('Started')

    if args.bathmet2:
        analyseBathMet2(args.bathmet2, args.reference)

    if args.bed:
        analyseBed(args.bed, args.reference)

    if args.bismark:
        analyseBismark(args.bismark, args.reference)

    if args.deepsignalplant:
        analyseDeepSignalPlant(args.deepsignalplant, args.reference)

