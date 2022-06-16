import argparse
import logging
from pathlib import Path
from time import asctime
from Bio import SeqIO
from tqdm import tqdm
import math
from collections import Counter
import pandas as pd
import pysam

parser = argparse.ArgumentParser(description='Process a fasta file for stats about the genome and C nucleotide. Add a GFF file for further analysis. You can add a nucleotide and the annotations you want.')
parser.add_argument('-f', '--fasta', type=Path, help='load fasta file')
parser.add_argument('-g', '--gff', type=Path, help='load gff file')
parser.add_argument('-n', '--nucleotide', default="C", type=str, help='give a nucleotide to get stats for')
parser.add_argument('-a', '--annotation', default="gene", type=str, help='give the annotations from the gff files you want to study')
parser.add_argument("-d", "--debug", dest="loglevel", help="For debugging and speed only. Use only a small part of the bam.", action="store_const", const=logging.DEBUG, default=logging.INFO)
args = parser.parse_args()

def getGenomeStats(fasta, outfile, nucleotide):
    # docu = http://biopython.org/DIST/docs/tutorial/Tutorial.html#sec53 chap 5
    logging.info('Starting getGenomeStats')
    total=0
    totalMb=0
    totalN=0
    with open(outfile, mode = "w") as out:
        for seq_ref in SeqIO.parse(fasta,"fasta"):
            print(seq_ref.id, file=out)
            logging.debug(repr(seq_ref.seq))
            print(len(seq_ref)/10**6,"Mb.", file=out)
            total+=len(seq_ref)
            totalMb+=len(seq_ref)/10**6
            compteN=0
            for nucl in tqdm((seq_ref)):
                if nucl == nucleotide:
                    compteN+=1
            totalN+=compteN
        print ("Le total est de",totalMb,"Mb.", file = out)
        print('Il y a ',totalN/total*100,'% de',nucleotide, file = out)
    logging.info('getGenomeStats finished')

def getNucleotidePosition (fasta,outfile, nucleotide):
    logging.info('Starting getNucleotidePosition')
    with open(outfile, mode = "w") as out:
        for seq_ref in SeqIO.parse(fasta,"fasta"):
            for index, nucl in tqdm(enumerate(seq_ref, start = 1)):
                if nucl == nucleotide:
                    print(seq_ref.id, index, file = out, sep="\t")
    logging.info('getNucleotidePosition finished')

def GetCPos(sequence):
    logging.info('Starting getNucleotidePosition')
    positions=[]
    for index, nucl in tqdm(enumerate(sequence, start = 1)):
        if nucl == 'C':
            positions.append(index)
            logging.debug(positions)
    logging.info('getGenomeStats finished')
    return positions
    
def GetCContext(sequence):
    logging.info('Starting getCContext')
    contextes=[]
    for index, nucl in tqdm(enumerate(sequence, start = 1)):
            if index+1 < len(sequence) :
                if nucl =='C':
                    contexte = sequence[index:index+2]
                    if contexte[0]=='G':
                         contexte= 'CG'
                    elif contexte[1]=='G':
                         contexte= 'CHG'
                    else : contexte= 'CHH'
                    contextes.append(contexte)
            else: pass
    logging.info('getCContext finished')
    return contextes

def dicoCompoGene(fasta):
    logging.info('Starting dicoCompoGene')
    genome = {}
    # genome[key] = valeur
    for seq_ref in SeqIO.parse(fasta,"fasta"):
        genome[seq_ref.id] = seq_ref.seq
    logging.info('dicoCompoGene finished')
    return genome

def getGff(gff, outfile, annotation) :
    logging.info('Starting getGff')
    with open (outfile, mode = "w") as out:
        with open(gff) as gff:
            for f in pysam.tabix_iterator(gff, parser = pysam.asGTF()):
                # Retrieve gene feature information

                if f.feature in annotation:
                    try:
                        # https://pysam.readthedocs.io/en/latest/api.html
                        attr_dict = dict(i.split("=") for i in f.attributes.split(";"))
                        name = attr_dict["Name"]
                        start = f.start
                        end = f.end
                        chromosome = f.contig
                        strand=f.strand

                    except ValueError:
                        raise Exception('Error. Cannot split the attributes in gff. Is the annotation file in GFF3 format?')
                    print(f.feature, name, chromosome,strand, start, end, sep="\t", file= out)
    logging.info('getGff finished')
    return f.feature, name, chromosome,strand, start, end

def getGeneComposition(fasta, gff, out, out2):
    logging.info('Starting getGeneComposition')
    genome=dicoCompoGene(fasta) 
    with open(gff, mode = "r") as gff:
        with open(out, 'w') as outfile:
        #print(gene.readlines())
            listetotal=[]
            nbrenucleotides={}
            print('Chr','Pos_Locale', 'Pos_Globale','ID','Taille','Contexte',sep='\t', file=outfile) #header
            for line in gff.readlines():
                if line.startswith('#'):
                    continue
                type, id, chr, strand, start, end = line.strip().split("\t")
                genesequence=genome[chr][int(start):int(end)]
                genesequencecont=genome[chr][int(start):int(end)+2]
                nbrenucleotides=Counter(genesequence) #fait un dico de la sequence
                list_C_relative = GetCPos(genesequence) #obtient la position des C relatifs aux elements
                list_contexte = GetCContext(genesequencecont)
                list_C_genomic = [int(index) + int(start) for index in list_C_relative]
                assert len(list_C_genomic) == len(list_C_relative), "Problem in list position" #permet de s'assurer de qqch, si faux message en ""
                taille=int(end)-int(start)+1
                for index_listC, index_list_G, contexte in zip(list_C_relative, list_C_genomic, list_contexte): #zip permet de lire plusieurs listes de mm taille en mm temps
                    print(chr, index_listC, index_list_G, id, taille, contexte, sep='\t', file=outfile) 
                liste=[type, id, chr, strand, start, end,nbrenucleotides['A'], nbrenucleotides['G'], nbrenucleotides['T'],nbrenucleotides['C'], taille, nbrenucleotides['C']/taille*100 ]
                listetotal.append(liste)
            tab=pd.DataFrame(listetotal, columns=['Annotation','ID','Chromosome', 'Brin','Start','End','Nbre_A', 'Nbre_G','Nbre_T','Nbre_C','Taille','Pourcentage_de_C'])
            logging.debug(tab)
            tab.sort_values(by=['Pourcentage_de_C'], ascending=False, inplace=True)
            tab.to_csv(out2, sep='\t', index=False)
    logging.info('getGeneComposition finished')
        

def getAnalyseRef(fasta=args.fasta, nucleotide=args.nucleotide, annotation=args.annotation, gff=args.gff) :
    getGenomeStats(fasta, nucleotide, 'Stats.txt')
    getNucleotidePosition(fasta, nucleotide+".txt",nucleotide)
    if not gff:
        logging.warning('Please load a gff file if you want a complete analysis.')
    else :
        getGff(gff, "GffTraite.txt", annotation)
        getGeneComposition(fasta,"GffTraite.txt" ,"RefPosC.txt", "Infos.txt")

def check_args(args):
    if args.fasta and not args.fasta.exists():
        logging.warning("Fasta input does not exist. Please give a valid file.")   
        exit()
    if args.gff and not args.gff.exists():
        logging.warning("Gff input does not exist. Please give a valid file.")
        exit()


if __name__ == '__main__':
    check_args(args)

    logging.basicConfig(filename='myapppos.log', filemode='w', format='%(asctime)s %(message)s', datefmt='%d/%m/%Y %H:%M:%S', level=args.loglevel)
    logging.info('Started')
    getAnalyseRef(fasta=args.fasta, nucleotide=args.nucleotide, annotation=args.annotation, gff=args.gff)
    logging.info('Finished')
 
    


    