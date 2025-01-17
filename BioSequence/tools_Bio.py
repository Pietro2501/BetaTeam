import argparse
import gzip
import os
import math
import sys

def GC_skew (seq: str) -> float:
    feature={}
    lunghezze = len(seq)
    for i in seq:
        for task in ('lun','perc_ACGT','GC'):
            if task=='lun':
                feature[i]= len(seq[i])
            elif task=='perc_ACGT':
                A=(seq[i].count('A')/lunghezze[i])*100
                C=(seq[i].count('C')/lunghezze[i])*100
                G=(seq[i].count('G')/lunghezze[i])*100
                T=(seq[i].count('T')/lunghezze[i])*100
                feature[i]= A,G,T,C
            elif task=='GC':
                G=seq[i].count('G')
                C=seq[i].count('C')
                skew= (G-C)/(G+C)*100
                feature[i]= skew
    print ('questo è il dizionario delle feature: ', feature)
    return feature

def check_file(file_path: str) -> bool:
    try:
        stat_file = os.stat(file_path)
        res = True
    except FileNotFoundError:
        sys.exit(f"{file_path} does not exist!")
    else:
        return res

def extract_info(filename):
    if filename.endswith('.gz') or filename.endswith('.gzip'):
        with gzip.open(filename,'rt') as file:
                contenuto = file.readlines()
                output_file = 'output_unzipped.fastq'
        with open(output_file,'w') as output:
            for line in contenuto:
                output.write(line)

    elif filename.endswith('.fastq') or filename.endswith('.fq'):
        with open(filename,'r') as file:
            c = file.readlines()
            output_file = 'output.fastq'
        with open(output_file,'w') as output:
            for line in c:
                output.write(line)
    else:
        print("Wrong File Format. Choose a FastQ or GZ file")

    return output_file


def generate_kmers(sequence, k):
    return {sequence[i:i+k] for i in range(len(sequence) - k + 1)}

def comp_rev (nameseq):
    nameseq = nameseq.replace('\n','')
    reverse= nameseq[::-1]
    table= str.maketrans('ACTG','TGAC')
    comp_rev= nameseq.translate(table)[::-1]
    return comp_rev

def median (val: list):
    lenght_val = len(val)
    if lenght_val%2 == 0:
        return (val[(lenght_val//2)-1]+val[lenght_val//2])/2
    else:
        return val[lenght_val//2]

def quartili(values: list):
    values.sort()
    sec_quart = median(values)
    values_cpy = values.copy()
    if len(values)%2 != 0:
        values_cpy.pop(len(values)//2)
    values_1_half = values_cpy[:len(values_cpy)//2]
    values_2_half = values_cpy[len(values_cpy)//2:]
    prim_quart = median(values_1_half)
    terz_quart = median(values_2_half)
    return prim_quart, sec_quart, terz_quart

def IQR(prim:float, terz:float):
    return terz - prim

def whiskers(prim: float, terz: float, IQR: float):
    lower = prim - 1.5*IQR
    higher = terz + 1.5*IQR
    return lower, higher

def calc_k_mer(k,list_seq:list):
    """
    :param k:
    :param contenuto:
    :return:
    """
    kmer_diz = {}
    for seq in list_seq:
        for i in range(0,len(seq)-k+1):
            kmer = seq[i:i+k]
            kmer_diz.setdefault(kmer, 0)
            kmer_diz[kmer] += 1
    return kmer_diz