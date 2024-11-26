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



 # Primer reverse: GGACTACNVGGGTWTCTAAT

#kozak_pattern2 = re.compile("(GCC)?([AG]CCATG[ACGT])")#GCC non molto conservata quindi metto in gruppo e la cerco 0 o più volte,opzionale
#for match in kozak_pattern2.finditer(dizionario['ENST00000003084.11']):
    #print(match.group(),match.span())

helico_pattern_fwd = regex.compile("GTGCCAGCMGCCGCGGTAA") #Primer fwd: GTGCCAGCMGCCGCGGTAA

def comp_rev (nameseq):
    nameseq = nameseq.replace('\n','')
    reverse= nameseq[::-1]
    table= str.maketrans('ACTG','TGAC')
    comp = nameseq.translate(table)
    comp_rev= nameseq.translate(table)[::-1]
    return comp,comp_rev

fw_comp =comp_rev('GTGCCAGCMGCCGCGGTAA')[0]
rev_comp = comp_rev('GGACTACNVGGGTWTCTAAT')[1]


# W -> A o T
# M -> A o C
# V -> A o C o G
# N -> A o T o C o G


fw_comp_pattern = regex.compile('CACGGTCG[TG]CGGCGCCATT')
rev_comp_pattern = regex.compile('ATTAGA[TA]ACCC[TGC][ACGT]GTAGTCC')



for match in fw_comp_pattern.finditer(helico['4762fdd895094939_1'],overlapped=True):
    print(match.group(), match.span())

for match in rev_comp_pattern.finditer(helico['4762fdd895094939_1'],overlapped=True):
    print(match.group(), match.span())



