__all__ = ["GetFastqPrefix","GetFastaPrefix"]

import re

def GetFastqPrefix(fastq):
    """
    Get the prefix of a fastq file
    :param fastq:
    :return:
    """
    return re.search(r'(.*?)_', fastq).group(1)

def GetFastaPrefix(fasta):
    """
    Get the prefix of a fasta file

    :param fasta:
    :return:
    """


    name = re.search(r'(\S+)\.fasta', fasta).group(1)

    name = str(name).split('_')

    return name[0]
