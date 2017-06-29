#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import numpy as np
import argparse
from Bio.Blast.Applications import NcbitblastnCommandline
from Bio.Blast import NCBIXML
from time import strftime

'''
GeCo (Gene Conservation) is a platform developed on python3 for finding gene
conservation along the refseq genomes of NCBI.

There is two ways to run the script depending if you have or not an output file
in xml of BLAST against Representative bacterial/archaeal genomes database of
NCBI ftp://ftp.ncbi.nlm.nih.gov/blast/db/Representative_Genomes.*.tar.gz

    1) If you have the output file of BLAST on xml format (from web service or
    locally) you have to run the script like this.

        python3 GeCo.py -g BLAST_OUTPUT.xml -b 75

    Where the flag "-b" is the threshold bit score for considerate a hit as a
    truly homologous.

    The output of GeCo is a raw tab-separated file and a matrix for a better
    analisis, and it contain Query, Hit and Bits data from BLAST file.

    gene_conservation_out contain a matrix with all hits as rows and all genes
    used as query, aditionaly this file has 2 columns more 'Gene conservation'
    and 'Bits Mean', that is the proportion of gene conservation for a single
    organism (Hit) and the mean of Bit score of all queries of this organism.

    2) If you do not have the output file of BLAST on xml format you must run
    BLAST, personally i recommend run BLAST like this:

        tblastn -query FASTA_AA.faa -db ~/BLASTDB/Representative_Genomes\
        -db_gencode 11 -matrix BLOSUM45 -num_threads 4  -evalue 0.001\
        -max_target_seqs 1000000 -max_hsps 1 -qcov_hsp_perc 75 -outfmt 5\
        -out tblastn_out.xml

    Then you have the input file for GeCo, so go to step 1

    2.a) Alternatively, GeCo can make a BLAST locally for you (but is better
    make the BLAST by yourself on terminal)

        python3 GeCo.py -blast FASTA_AA.faa -db ~/BLASTDB/Representative_Genomes\
        -o tblastn_out.xml

    Once you have the "tblastn_out.xml" you must to go to step 1

Dependencies:
    
    pip3 install pandas
    pip3 install numpy
    pip3 install biopython
'''

__author__ = 'Enzo Guerrero-Araya (biologoenzo@gmail.com)'
__version__ = '0.1'
__date__ = 'June 14, 2016'


TIME = strftime('%H%M%S-%d%m%y')


def make_blastp(QUERY, DB, outfile):
    tblastn = NcbitblastnCommandline(query=QUERY, db=DB, outfmt=5,
                                     out=outfile, db_gencode=11,
                                     matrix='BLOSUM45',
                                     evalue=0.001, qcov_hsp_perc=75,
                                     num_alignments=1000000, num_threads=4,
                                     max_hsps=1)
    stdout_p, stderr_p = tblastn()
    print('BLAST DONE [{} file was created on this folder]'.format(outfile))


def _get_mean(names, table):
    x = [table[name] for name in names]
    return sum(x) / len(x)


def _to_bin(x):
    if x > 0:
        return 1
    else:
        return 0


def _get_percentage(names, table):
    x = []
    table_cp = table.copy()
    for name in names:
        table_cp[name] = list(table_cp[name].get_values())
        table_cp[name] = list(map(_to_bin, table_cp[name]))
        x.append(table_cp[name])
    return sum(x) / len(names)


def get_genome_name(blastxml, BITS=75):
    blast_result = NCBIXML.parse(blastxml)
    tsv = open('rawdata_out_{}.tsv'.format(TIME), 'w')
    print('Query', 'Hit', 'Bits', sep='\t', file=tsv)
    for i, blast_record in enumerate(blast_result):
        query = blast_record.query
        for aln in blast_record.alignments:
            hit = aln.title
            bits = aln.hsps[0].bits
            if bits >= BITS:
                print(query, hit, bits, sep='\t', file=tsv)
    tsv.close()
    table = pd.read_table('rawdata_out_{}.tsv'.format(
        TIME), engine='python', sep='\t', header=0)
    table = table.pivot_table(
        values='Bits', index='Hit', columns='Query', fill_value=0)
    table['Gene conservation'] = _get_percentage(
        list(table.columns.values), table)
    table['Bits Mean'] = _get_mean(list(table.columns.values)[:-1], table)
    table = table.sort_values(
        ['Gene conservation', 'Bits Mean'], ascending=False)
    table.to_csv('gene_conservation_out_{}.tsv'.format(TIME), sep='\t')
    print('GeCo DONE [gene_conservation_out_{}.tsv file was created on this folder]'.format(TIME))


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Make a GeCo for finding the \
                        conservation of a gene set along the refseq genomes \
                        of NCBI')
    parser.add_argument('-g', '--geco', dest='blastxml',
                        help='Output file of BLAST on xml format',
                        type=argparse.FileType('r'), metavar='BLASTXML')
    parser.add_argument('-b', '--bits', dest='bits',
                        help='The bits threshold for considerate a hit as a \
                        truly homologous [be careful with false positives] \
                        (default: 75)', default=75)
    parser.add_argument('-blast', dest='fastas',
                        help='Fasta file of aminoacid sequences to \
                        use as query',
                        metavar='FASTA_AA')
    parser.add_argument('-db', dest='db',
                        help='BLAST DB of Representative Genomes of NCBI')
    parser.add_argument('-o', '--out-file', dest='outfile',
                        help='File name for BLAST output')

    args = parser.parse_args()
    if (args.outfile is None or args.db is None) and args.fastas is not None:
        parser.error("-o AND -db is requires for make a BLAST")
    else:
        db = args.db
        outfile = args.outfile
        fastas = args.fastas

    if args.blastxml is not None:
        blastxml = args.blastxml
        bits = args.bits
        get_genome_name(blastxml, bits)
    elif args.fastas is not None:
        make_blastp(fastas, db, outfile)
    else:
        parser.print_help()
