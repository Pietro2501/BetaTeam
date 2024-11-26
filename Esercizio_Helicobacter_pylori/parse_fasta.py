import os

def parse_fasta (file_name):
    os.path.exists('file_name.fa')
    with open('file_name.fa', 'r') as f:
        contenuto = f.readlines()

    seq={}
    for i in contenuto:
        if i.startswith('>'):
            name=i.strip().split()[0][1:]
            seq[name]=''
        else:
            seq[name]=seq[name]+i.strip()
    return seq    




