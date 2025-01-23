from BioSequence import tools_Bio as t
import pickle
import gzip

def parse_fastq(file_obj: str, offset=33) -> dict:
    """
    Parses a FASTQ file and returns a dictionary representation of the sequences,
    quality scores and ASCII quality scores.

    :param file_obj: A file object containing the FASTQ data to be parsed.
    :type file_obj: str
    :param offset: Offset value used for ASCII quality score conversion.
                   Should only be either 33 or 64.
    :type offset: int
    :return: A dictionary where the keys are sequence identifiers (headers starting
             with '@') and values are sub-dictionaries containing sequence strings,
             ASCII quality scores, and numeric quality scores.
    :rtype: dict or None
    """
    with open(file_obj, "r") as file_obj:
        def check_offset(offset: int) -> bool:
            res = False
            if type(offset) == int and offset in [33, 64]:
                res = True
            return res

        def check_coherence(acc: str, plus: str) -> bool:
            """
            Validates the coherence of a FASTQ record by checking the format of
            the accession line and the "+" line.

            :param acc: The accession line from the FASTQ record.
            :type acc: str
            :param plus: The "+" line from the FASTQ record.
            :type plus: str
            :return: True if both the accession and the "+" lines are in the
                correct format, otherwise False.
            :rtype: bool
            """
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
                seq, plus, qual = file_obj.readline().strip(), file_obj.readline().strip(), file_obj.readline().strip()
                if not check_coherence(acc, plus):
                    print("Il file è corrotto!!!")
                    fastq_dict = None
                    break
                if seq and qual:
                    fastq_dict[acc]["seq"] = seq
                    fastq_dict[acc]["ASCII_qual"] = qual
                    fastq_dict[acc]["qual"] = [ord(q) - offset for q in fastq_dict[acc]["ASCII_qual"]]
                else:
                    print("Il file è corrotto!!!")
                    fastq_dict = None
                    break
                line = file_obj.readline()
            return fastq_dict
        elif not check_offset(offset):
            print(f"l'offset indicato non è numerico o è un valore differente da 33 o 64!!!")


def exp_err(qual, is_ascii=False):
    """
    Calculate the sum of expected error probabilities from quality scores.

    This function computes the expected error based on provided quality scores.
    The quality scores can either be in ASCII-encoded form (default) or as raw
    numerical values.

    :param qual: The quality scores. If `is_ascii` is True, they are
        expected to be ASCII-encoded characters. Otherwise, they should
        be numerical values.
    :type qual: list or str
    :param is_ascii: Whether the input `qual` is ASCII-encoded. Defaults
        to False.
    :type is_ascii: bool
    :return: The sum of expected error probabilities.
    :rtype: float
    """
    if is_ascii:
        ee = sum(map(lambda q: 10 ** (q / -10), map(lambda x: ord(x) - 33, qual)))
    else:
        ee = sum(map(lambda q: 10 ** (q / -10), qual))
    return ee


def dinamic_trimming(dict_fastq: dict, threshold=25, window=15):
    """
    Performs dynamic trimming of sequencing reads based on quality scores.

    This function takes a dictionary containing sequencing reads and their quality
    information. It performs a dynamic trimming process by evaluating the average
    quality score of a sliding window. When the average quality falls below a given
    threshold, the function trims the read from the point where quality criteria
    are no longer met.

    :param dict_fastq: Dictionary containing sequencing reads and their associated
        quality scores.
    :param threshold: Integer threshold below which the average quality score is
        considered too low. Default is 25.
    :param window: Integer size of the sliding window used to compute the average
        quality score. Default is 15.
    :return: The updated dictionary after processing.
    :rtype: dict
    """
    for key in dict_fastq.keys():
        quality = dict_fastq[key]["qual"]
        list_mean = quality[:window]
        if sum(list_mean)/len(list_mean) < threshold:
            for value in range(len(list_mean), len(quality)):
                if (sum(list_mean)+quality[value])/(len(list_mean)+1) < threshold:
                    list_mean.append(quality[value])
                else:
                    break
            dict_fastq[key].update({"seq": dict_fastq[key]["seq"][len(list_mean):]})
            dict_fastq[key].update({"qual": dict_fastq[key]["qual"][len(list_mean):]})
            dict_fastq[key].update({"ASCII_qual": dict_fastq[key]["ASCII_qual"][len(list_mean):]})
    return dict_fastq

def dict_filter(_dict_fastq: dict, k: int):
    """
    The function iterates through the dictionary, calculates an expected error value
    for each entry's quality scores, and deletes entries where the error value
    exceeds the given threshold.

    :param _dict_fastq: Dictionary containing FASTQ data where the key is typically
        the identifier of a read and the value is a dictionary containing read
        attributes.
    :type _dict_fastq: dict
    :param k: The error threshold. Entries in `_dict_fastq` with expected error
        values greater than this threshold will be removed.
    :type k: int
    :return: A dictionary with filtered FASTQ entries where all remaining entries meet
        the error threshold defined by `k`.
    :rtype: dict
    """
    list_keys = list(_dict_fastq.keys())
    for key in list_keys:
        expect_err = exp_err(_dict_fastq[key]["qual"])
        if expect_err > k:
            del _dict_fastq[key]
    return _dict_fastq

def rev_comp(seqlist: list):
    """
    Generate the reverse complement sequences for a given list of DNA sequences, append them to the
    original list, and return the combined list.

    :param seqlist: A list of DNA sequences for which reverse complements are to be calculated
    :type seqlist: list
    :return: A combined list of original sequences and their reverse complements
    :rtype: list
    """
    table= str.maketrans('ACTG','TGAC')
    total_list = [seq.translate(table)[::-1] for seq in seqlist] + seqlist
    return total_list

if __name__ == "__main__":
    file_f= t.extract_info(r"Phoecicola_V_sim1.fq.gz")
    dict_fastq_f= parse_fastq(file_f)
    # print (f"il file R1 contiene {len(dict_fastq_f)} sequenze")
    # dict_filtered_f = dinamic_trimming(dict_fastq_f)
    # dict_double_filtered_f = dict_filter(dict_filtered_f, 3)
    # seq_filtered_f = [dict_double_filtered_f[key]["seq"] for key in dict_double_filtered_f.keys() if len(dict_double_filtered_f[key]["seq"]) > 75]
    seq_f = [dict_fastq_f[key]["seq"] for key in dict_fastq_f.keys()]
    print (f"andremo a lavorare con {len(seq_f)} sequenze del file R1")

    file_r= t.extract_info(r"Phoecicola_V_sim2.fq.gz")
    dict_fastq_r= parse_fastq(file_r)
    # print (f"il file R2 contiene {len(dict_fastq_r)} sequenze")
    # dict_filtered_r = dinamic_trimming(dict_fastq_r)
    # dict_double_filtered_r = dict_filter(dict_filtered_r, 3)
    # seq_filtered_r = [dict_double_filtered_r[key]["seq"] for key in dict_double_filtered_r.keys() if len(dict_double_filtered_r[key]["seq"]) > 75]
    seq_r = [dict_fastq_r[key]["seq"] for key in dict_fastq_r.keys()]
    print (f"andremo a lavorare con {len(seq_r)} sequenze del file R2")

    mucchio_selvaggio = seq_f + seq_r
    delirio_totale = rev_comp(mucchio_selvaggio)
    print(f"il grafo verrà costruito utilizzando {len(delirio_totale)} sequenze")


    with gzip.open("kmer_17_NEW_PhoeVul_num1filtered.pkl.gz", "wb") as f:
        pickle.dump(t.calc_k_mer(17, delirio_totale), f)

