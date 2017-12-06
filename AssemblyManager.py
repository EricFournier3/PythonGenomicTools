
import os
import sys
import re
import argparse
import yaml
from Bio import SeqIO

'''
######################################################################
# This script is used to perform filtration or compute statistics on #
# a set of Spades assembly files.                                             #
######################################################################

Command line example for filtration:
-----------------------------------

python AssemblyManager.py --assembdir '/home/eric/ProjetPersonnel/Programmation/ProjetPython/NGS/TEST/' --do 'filter' --yaml '/home/eric/ProjetPersonnel/Programmation/ProjetPython/NGS/TEST/myyaml.yaml' --length 1000 --cov 5

'''


#Command line options
parser = argparse.ArgumentParser(prog="Assembly Manager")

#Base directory with Spades fasta assembly
parser.add_argument('--assembdir',metavar='[Required : Path to fasta assembly]',required=True,nargs=1)

#User can choose to filter assembly files or compute statistics
parser.add_argument('--do',nargs=1,type=str,choices=['filter','stat'],metavar='[Required : Type of treatment]',required=True)

#Path to the yaml file with specimen list
parser.add_argument('--yaml',type=str,metavar='[Required : Yaml file with specimen id]',nargs=1,required=True)

#Contig length threshold
parser.add_argument('--length',type=int,metavar='[Required : length threshold for filtration]',required=True,nargs=1)

#Contig coverage threshold
parser.add_argument('--cov',type=int,metavar='[Required : coverage threshold for filtration]',required=True,nargs=1)

args=parser.parse_args()


class Filter():
    """
    Class for Spades assemblies filtration
    """

    def __init__(self):

        #Path to yaml file with specimen id list
        self.yaml_file_path = args.yaml[0]

        #Base directory with assembly fasta file
        self.assembly_dir = args.assembdir[0]

        if (not self.assembly_dir.endswith('/')):
            self.assembly_dir = self.assembly_dir + '/'

        #Output directory for filtrated assembly files
        self.out = self.assembly_dir + r"Filtred/"

        #Create the output directory
        try:
            os.mkdir(self.out)
        except Exception as e:
            print "In __init__"
            print e
            exit(0)

        #Contig Length and coverage threshold
        self.length_thresh = args.length[0]
        self.cov_thresh = args.cov[0]

    def ParseYamlFile(self):
        """
        Build the specimen list
        :return:
        """

        #yaml file with specimen list
        self.yam_file = open(self.yaml_file_path).read()

        #number of specimen
        self.nb_spec = yaml.load(self.yam_file)[0]

        #Specimen list
        self.SpecList = yaml.load(self.yam_file)[1:self.nb_spec + 1]

    def DoFiltration(self,specimen):
        """
        Perform filtration of a Spades assembly files
        :param specimen:
        :return:
        """

        #Contigs to keep according to the coverage and length threshold
        rec_to_keep = []

        try:

            for my_rec in SeqIO.parse("{0}{1}{2}".format(self.assembly_dir,specimen,'.fasta'), "fasta"):

                id_parse = re.search(r'NODE_\d+_length_(\d+)_cov_(\S+)_ID_',my_rec.id)
                contig_length = id_parse.group(1)
                contig_cov = id_parse.group(2)

                if((int(contig_length) >= self.length_thresh)  and (float(contig_cov) >= self.cov_thresh)):
                    rec_to_keep.append(my_rec)

            SeqIO.write(rec_to_keep,self.out + specimen + "_Over_" + str(self.length_thresh) + "pb_" + str(self.cov_thresh) + "X"  + '.fasta','fasta')

        except Exception as e:
            print "In DoFiltration"
            print (e)


class AssembStat():

    def __init__(self):
        pass
        #TODO to complete 



if (args.do[0] == 'filter'):
    my_filter = Filter()
    my_filter.ParseYamlFile()
    map(my_filter.DoFiltration,my_filter.SpecList)
else:
    my_stat = AssembStat()

