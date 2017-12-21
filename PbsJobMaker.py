import sys
import argparse
import os
import re


'''
##################################################################################################################
# This script produce one Briaree server pbs job file for every paired end fastq.gz found in the input directory #
##################################################################################################################

Command example:
- - - - - - - -
python PbsJobMaker.py  --fastqdir '/home/eric/ProjetPersonnel/Programmation/ProjetPython/NGS/TEST/FASTQ'  
--outdir '/home/eric/ProjetPersonnel/Programmation/ProjetPython/NGS/TEST/OUTPBS' --wt 2 -n 1 --ppn 12 -q 'normale'  
--gtp '/home/eric/ProjetPersonnel/Programmation/ProjetPython/NGS/InGit/GenomicTools/'  --cqfastqdir '20171201_VNO' 
--cquser 'fournie1' -k 21 33 55 71 111 --cqspades 'SPAdes-3.11.1-Linux'

'''


#Command line options
parser = argparse.ArgumentParser(prog="Pbs Job Maker")

#Local Base directory with fastq.gz files
parser.add_argument('--fastqdir',type=str,metavar='[Required : Path to local fastq files]',required=True,nargs=1)

#Calcul Quebec fastq directory
parser.add_argument('--cqfastqdir', type=str, metavar='[Required : Calcul Quebec fastq directory on Briaree scratch]',required=True,nargs=1)

#Calcul Quebec user name
parser.add_argument('--cquser', type=str, metavar='[Required : Calcul Quebec user name on Briaree server]',required=True,nargs=1)

#Spades version
parser.add_argument('--cqspades', type=str, metavar='[Required : Calcul Quebec spades version on Briaree scratch]',required=True,nargs=1)

#K-mer list for spades assembly
parser.add_argument('-k', type=str,metavar='[Required : Kmer list for spades assembly]', required=True,nargs='+')

#Pbs output directory
parser.add_argument('--outdir',type=str,metavar='[Required : Path to output directory]',required=True,nargs=1)

#Job Wall time
parser.add_argument('--wt', type=int, metavar='[Required : Wall time]',required=True,nargs=1)

#Number of nodes
parser.add_argument('-n',type=int,metavar='[Required : number of nodes]',required=True,nargs=1)

#Number of core per node
parser.add_argument('--ppn',type=int,metavar='[Required : processors per node]',required=True,nargs=1)

#Queue name
parser.add_argument('-q',type=str,metavar='[Required : queue name]',choices=['normale','test','courte'],required=True,nargs=1)

#Path to the GenomicTools modules
parser.add_argument('--gtp',type=str,metavar='[Required : path to GenomicTools modules]',required=True,nargs=1)

args=parser.parse_args()

if (not os.path.exists(args.gtp[0])):
    print "No GenomicsTools modules"
    exit(0)

sys.path.append(args.gtp[0])


import PythonUtils
from LinuxUtils import Bash


class PbsWriter(object):
    """
    Class to build Pbs file
    """
    def __init__(self):

        #Init with command line arguments
        self.fastq_dir = args.fastqdir[0]
        self.out_dir = args.outdir[0]
        self.walltime = args.wt[0]
        self.node = args.n[0]
        self.ppn = args.ppn[0]
        self.queue = args.q[0]
        self.cquser = args.cquser[0]
        self.kmer = args.k

        #Spades bin localization in Briaree $SCRATCH
        self.cqspades = '/RQexec/{0}/{1}/bin'.format(self.cquser,args.cqspades[0])

        #Path to fastq files in Briaree $SCRATCH
        self.cqfastq_dir = '/RQexec/{0}/{1}/'.format(self.cquser,str(args.cqfastqdir[0]).strip('/'))

        #Path where to copy fasta output in Briaree $SCRATCH
        self.cqfasta_out = self.cqfastq_dir + "Assemblies/"

        #Spades output path
        self.cqspades_out = self.cqfastq_dir + "SPADES/"


        if (not str(self.out_dir).endswith('/')):
            self.out_dir += '/'

        #List of fastq files
        self.fastq_list = []


    def GetHeader(self,specimen):
        """
        Build the Pbs header for this specimen
        :param specimen:
        :return:
        """

        s = "#!/bin/sh\n"
        s += "#\n"
        s += "#PBS -l walltime={0}:00:00\n".format(self.walltime)
        s += "#PBS -l nodes={0}:ppn={1}\n".format(self.node, self.ppn)
        s += "#PBS -o {0}.dat\n".format(specimen)
        s += "#PBS -j oe\n"
        s += "#PBS -W umask=022\n"
        s += "#PBS -r n\n"
        s += "#PBS -q {0}\n\n\n".format(self.queue)

        return s


    def GetFastqList(self):
        """
        Get a fastq.gz file name list from the local fastq directory
        :return:
        """

        self.fastq_list = Bash.GetFastqInDir(self.fastq_dir,'gz')


    def CreatePbsDir(self):
        """
        Make the Pbs output directory
        :return:
        """

        if (os.path.exists(self.out_dir)):
            pass
        else:
            os.mkdir(self.out_dir)

    def GetCore(self,fastq):
        """
        Build the Pbs bash core for this specimen
        :param fastq:
        :return:
        """

        s = "FASTADIR=\"{0}\"\n".format(self.cqfasta_out)
        s += "FASTQDIR=\"{0}\"\n".format(self.cqfastq_dir)
        s += "SPADESDIR=\"{0}{1}/\"\n".format(self.cqspades_out,PythonUtils.GetFastqPrefix(fastq))
        s += "p1_1=\"{0}\"\n".format(fastq)
        s += "p1_2=\"{0}\"\n\n".format(re.sub('_R1_','_R2_',fastq))
        s += "cd {0}\n".format(self.cqspades)
        s += "./spades.py -k {0} --pe1 {1} --pe2 {2} --careful -o {3}\n\n".format(",".join(self.kmer),"$FASTQDIR$p1_1","$FASTQDIR$p1_2","$SPADESDIR")
        s += "cp $SPADESDIR\"contigs.fasta\" $FASTADIR\"{0}.fasta\"\n\n".format(PythonUtils.GetFastqPrefix(fastq))

        return s

    def WritePbsFile(self,fastq):
        """
        Write a Pbs file for this specimen
        :param fastq:
        :return:
        """

        if(str(fastq).find('_R1_') > -1):

            #Extract the specimen name from the R1 fastq file name
            specimen = PythonUtils.GetFastqPrefix(fastq)

            #The Pbs output file
            handler = open(self.out_dir + specimen + '.pbs','w')
            handler.write(self.GetHeader(specimen))
            handler.write(self.GetCore(fastq))

            handler.close()


my_pbs_writer = PbsWriter()
my_pbs_writer.CreatePbsDir()
my_pbs_writer.GetFastqList()


#Build a Pbs file for every specimen
map(my_pbs_writer.WritePbsFile,my_pbs_writer.fastq_list)

