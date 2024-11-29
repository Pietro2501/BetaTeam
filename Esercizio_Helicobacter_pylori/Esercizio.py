import os
import regex

def parse_fasta (file_name):
    os.path.exists(file_name)
    with open(file_name, 'r') as f:
        contenuto = f.readlines()

    seq={}
    for i in contenuto:
        if i.startswith('>'):
            name=i.strip().split()[0][1:]
            seq[name]=''
        else:
            seq[name]=seq[name]+i.strip()
    return seq

helico = parse_fasta('Helicobacter_pylori_ATCC_700392.fa')
#print(helico)

def fn_comp_rev (nameseq):
    nameseq = nameseq.replace('\n','')
    reverse= nameseq[::-1]
    table= str.maketrans('ACTG','TGAC')
    comp = nameseq.translate(table)
    comp_rev= nameseq.translate(table)[::-1]
    return comp,comp_rev

fw = 'GTGCCAGCMGCCGCGGTAA'
rev_comp = fn_comp_rev('GGACTACNVGGGTWTCTAAT')[1]
# print(rev)

# W -> A o T
# M -> A o C
# V -> A o C o G
# N -> A o T o C o G


fw_comp_pattern = regex.compile('GTGCCAGC[AC]GCCGCGGTAA')
rev_comp_pattern = regex.compile('ATTAGA[AT]ACCC[TCG][ATCG]GTAGTCC')

list_fwd=[]
list_rev=[]

for match in fw_comp_pattern.finditer(helico['4762fdd895094939_1'],overlapped=True):
    # print(match.group(), match.span())
    list_fwd.append(match.span()[0])

for match in rev_comp_pattern.finditer(helico['4762fdd895094939_1'],overlapped=True):
    # print(match.group(), match.span())
    list_rev.append(match.span()[1])

for j in list_rev:
    for i in list_fwd:
        if j-i >= 250 and j-i <= 350:
            print(f"I possibili ampliconi associati ai primer usati vanno da {i} a {j} e la loro lunghezza è pari a: {j-i}")

