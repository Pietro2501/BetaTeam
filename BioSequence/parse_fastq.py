import sys
import os

def parse_fastq(file_obj, offset=33) -> dict:
    """

    :param file_obj:
    :param offset:
    :return:
    """
    def check_offset(offset: int) -> bool:
        res = False
        if type(offset) == int and offset in [33,64]:
            res = True
        return res
    def check_coherence(acc: str, plus: str) -> bool:
        res = False
        if acc.startswith("@") and plus.startswith("+"):
            res = True
        return res
    if check_offset(offset):
        fastq_dict = {}
        line = file_obj.readline()
        while line:
            acc = line.strip()
            fastq_dict.setdefault(acc, {})
            # qui andrebbe inserito un try/expcet per verificare che effettivamente tutte le linee siano lette
            seq, plus, qual = file_obj.readline().strip(), file_obj.readline().strip(), file_obj.readline().strip()
            if not check_coherence(acc, plus):
                print("Il file è corrotto!!!")
                fastq_dict = None
                break
            if seq and qual:
                fastq_dict[acc]["seq"] = seq
                fastq_dict[acc]["ASCII_qual"] = qual
                fastq_dict[acc]["qual"] = [ord(q) - offset for q in fastq_dict[acc]["ASCII_qual"]]
                #fastq_dict[acc]["Pe"] = [10 ** (q / -10) for q in fastq_dict[acc]["qual"]]
            else:
                print("Il file è corrotto!!!")
                fastq_dict = None
                break
            line = file_obj.readline()
        return fastq_dict
    elif not check_offset(offset):
        print(f"l'offset indicato non è numerico o è un valore differente da 33 o 64!!!")

def exp_err(qual, is_ascii=False):
    if is_ascii:
        ee = sum(map(lambda q: 10** (q/-10), map(lambda x: ord(x) -33 ,qual)))
    else:
        ee = sum(map(lambda q: 10** (q/-10), qual))
    return ee