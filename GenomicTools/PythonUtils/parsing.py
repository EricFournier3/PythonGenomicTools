__all__ = ["GetFastqPrefix"]

import re

def GetFastqPrefix(fastq):
    """
    Get the prefix of a fastq file
    :param fastq:
    :return:
    """
    return re.search(r'(.*?)_', fastq).group(1)
