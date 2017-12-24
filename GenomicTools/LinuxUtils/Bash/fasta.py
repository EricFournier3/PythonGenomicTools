__all__ = ["GetFastaInDir"]


import subprocess
import os

def GetFastaInDir(path):
    """
    Return a list of fasta files in this path
    :param path:
    :param ext:
    :return:
    """

    if(not str(path).endswith('/')):
        path += '/'

    if(os.path.isdir(path)):

        command = "basename -a $(ls -1 {0}*.{1})".format(path,'fasta')

        proc = subprocess.Popen([command],stdout=subprocess.PIPE, shell=True)

        (out, err) = proc.communicate()

        fasta_list = out.split('\n')
        return sorted(fasta_list[0:len(fasta_list) - 1])

    else:
        print "Directory doesn't exist"
        exit(0)