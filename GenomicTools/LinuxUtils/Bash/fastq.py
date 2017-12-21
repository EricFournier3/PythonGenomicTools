__all__ = ["GetFastqInDir"]


import subprocess
import os


def GetFastqInDir(path,ext):
    """
    Return a list of fastq files having extension ext in path
    :param path:
    :param ext:
    :return:
    """

    if(not str(path).endswith('/')):
        path += '/'

    if(os.path.isdir(path)):

        command = "basename -a $(ls -1 {0}*.{1})".format(path,ext)

        proc = subprocess.Popen([command],stdout=subprocess.PIPE, shell=True)

        (out, err) = proc.communicate()

        fastq_list = out.split('\n')
        return sorted(fastq_list[0:len(fastq_list) - 1])

    else:
        print "Directory doesn't exist"
        exit(0)



