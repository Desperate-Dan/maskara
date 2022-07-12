#!/usr/bin/env python
from Bio import SeqIO
import itertools
import os
import argparse
import pysam
import sys

#Main body of the function

def runner(args):
    input_abspath = os.path.abspath(args.input_file)
    print(input_abspath)

def main():
    parser = argparse.ArgumentParser()
    '''
    Parse the command line arguments.
    '''
    parser = argparse.ArgumentParser(description='Creates a coverage mask to apply to your lovely consensus fasta.')

    required_group = parser.add_argument_group('Required')
    required_group.add_argument('-f', '--file', dest='input_file', required=True,
                            help='Path to the BAM/SAM file you want to create a mask for.')

    optional_group = parser.add_argument_group('Optional')
    optional_group.add_argument('-d', '--depth', dest='depth', default="20",
                            help='If coverage is below this it will be masked.')
       
    args = parser.parse_args()
    runner(args)    
    
    
    
    
    
if __name__ == "__main__":
    main()
