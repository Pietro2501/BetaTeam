import os
from idlelib.help import HelpWindow


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
print(helico)