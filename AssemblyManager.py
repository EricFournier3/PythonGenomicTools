
import os
import sys
import re
import argparse

from Bio import SeqIO

'''
######################################################################
# This script is used to perform filtration or compute statistics on #
# a set of Spades assembly files.                                    #
######################################################################

Command line example for filtration:
-----------------------------------

python AssemblyManager.py --assembdir '/home/eric/ProjetPersonnel/Programmation/ProjetPython/NGS/TEST/' --do 'filter' 
--gtp '/home/eric/ProjetPersonnel/Programmation/ProjetPython/NGS/InGit/GenomicTools/' --length 1000 --cov 5

'''


#Command line options
parser = argparse.ArgumentParser(prog="Assembly Manager")

#Base directory with Spades fasta assembly
parser.add_argument('--assembdir',metavar='[Required : Path to fasta assembly]',required=True,nargs=1)

#User can choose to filter assembly files or compute statistics
parser.add_argument('--do',nargs=1,type=str,choices=['filter','stat'],metavar='[Required : Type of treatment]',required=True)


#Contig length threshold
parser.add_argument('--length',type=int,metavar='[Required : length threshold for filtration]',required=True,nargs=1)

#Contig coverage threshold
parser.add_argument('--cov',type=int,metavar='[Required : coverage threshold for filtration]',required=True,nargs=1)

#Path to the GenomicTools modules
parser.add_argument('--gtp',type=str,metavar='[Required : path to GenomicTools modules]',required=True,nargs=1)

args=parser.parse_args()

if (not os.path.exists(args.gtp[0])):
    print "No GenomicsTools modules"
    exit(0)

sys.path.append(args.gtp[0])

import PythonUtils
from LinuxUtils import Bash


class AssemblyCommon(object):
    """
    Base class for class Filter and class AssembStat
    """

    def __init__(self):


        # Base directory with assembly fasta file
        self.assembly_dir = args.assembdir[0]

        # Contig Length and coverage threshold
        self.length_thresh = args.length[0]
        self.cov_thresh = args.cov[0]

        if (not self.assembly_dir.endswith('/')):
            self.assembly_dir = self.assembly_dir + '/'


    def MakeSpecList(self):

        fasta_list  = Bash.GetFastaInDir(self.assembly_dir)

        self.SpecList = [PythonUtils.GetFastaPrefix(x) for x in list(fasta_list) ]



class Filter(AssemblyCommon):
    """
    Class for Spades assemblies filtration
    """

    def __init__(self):

        #Base class
        AssemblyCommon.__init__(self)

        #Output directory for filtrated assembly files
        self.out = self.assembly_dir + r"Filtred"+ "_Over_" + str(self.length_thresh) + "pb_" + str(self.cov_thresh) + "X"  + "/"

        #Create the output directory
        try:
            os.mkdir(self.out)
        except Exception as e:
            print "In __init__"
            print e
            exit(0)


    def DoFiltration(self,specimen):
        """
        Perform filtration of a Spades assembly files
        :param specimen:
        :return:
        """

        #Contigs to keep according to the coverage and length threshold
        rec_to_keep = []

        try:
            #For each record in the assembly file
            for my_rec in SeqIO.parse("{0}{1}{2}".format(self.assembly_dir,specimen,'.fasta'), "fasta"):

                #Fasta header parsing
                id_parse = re.search(r'NODE_\d+_length_(\d+)_cov_(\S+)_ID_',my_rec.id)
                contig_length = id_parse.group(1)
                contig_cov = id_parse.group(2)

                #Keep the record only if his length and coverage are higher than threshold
                if((int(contig_length) >= self.length_thresh)  and (float(contig_cov) >= self.cov_thresh)):
                    rec_to_keep.append(my_rec)

            #Save the new filtered assembly
            SeqIO.write(rec_to_keep, self.out + specimen + '.fasta', 'fasta')

        except Exception as e:
            print "In DoFiltration"
            print (e)


class AssembStat(AssemblyCommon):
    """
    Compute assemblies statistics
    """

    def __init__(self):

        #Base class
        AssemblyCommon.__init__(self)

        # Output directory for statistics files
        self.out = self.assembly_dir + r"Stat/"

        # Create the output directory
        try:
            os.mkdir(self.out)
        except Exception as e:
            print "In __init__"
            print e
            exit(0)


    def CompileLengthAndCov(self,specimen):
        """
        Compilation of contig length and coverage
        for this specimen
        :param specimen:
        :return:
        """

        #Mean GC percent
        mean_gc = lambda gc,length : round(gc * 100.0 / length, 0)

        #Mean weigthed coverage
        mean_weighted_coverage =  lambda covlenList,length : round(sum([float(x)/float(length) for x in covlenList]),0)

        #Output statistic file
        stat_file = self.out + "AssembStat_" + specimen + ".txt"

        #Total contig length
        total_length = 0

        #Total GC content
        total_gc = 0

        #List of coverage x contig length
        covlen_list = []

        #Write the compilation in the output file
        with open(stat_file, 'w') as writef:

            #Header
            writef.write("NODE\tCoverage\tLength\n")

            try:

                #Parse every record
                for my_rec in SeqIO.parse("{0}{1}{2}".format(self.assembly_dir, specimen, '.fasta'), "fasta"):

                    #Fasta header parsing
                    id_parse = re.search(r'NODE_(\d+)_length_(\d+)_cov_(\S+)_ID_', my_rec.id)

                    #Node ID
                    node = id_parse.group(1)

                    #Contig length
                    contig_length = id_parse.group(2)

                    #Contig coverage
                    contig_cov = id_parse.group(3)

                    #Keep this contig for compilation only if his length and his coverage is higher than threshold
                    if(int(contig_length) >= self.length_thresh and  int(re.search(r'(\d+)\.?',contig_cov).group(1)) >= self.cov_thresh):
                        writef.write(node + "\t" + contig_cov + "\t" + contig_length + "\n")

                        total_length += int(contig_length)
                        total_gc += sum(my_rec.seq.count(x) for x in ['G', 'C', 'g', 'c', 'S', 's'])
                        covlen_list.append(int(re.search(r'(\d+)\.?',contig_cov).group(1))*int(contig_length))

                writef.write("\n")

                #Global Statistics
                writef.write("Total length\tPercent GC\tMean coverage\n")
                writef.write(str(total_length) + '\t' + str(mean_gc(total_gc,total_length)) + '\t' + str(mean_weighted_coverage(covlen_list,total_length)))


            except Exception as e:
                print "In CompileLengthAndCov"
                print (e)

        writef.close()


#Contigs filtration
if (args.do[0] == 'filter'):
    my_filter = Filter()
    my_filter.MakeSpecList()
    map(my_filter.DoFiltration,my_filter.SpecList)
#Statistics calculation
else:
    my_stat = AssembStat()
    my_stat.MakeSpecList()
    map(my_stat.CompileLengthAndCov, my_stat.SpecList)

