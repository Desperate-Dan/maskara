#!/usr/bin/env python3

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import itertools
import os
import argparse
import pysam
import sys

#Main body of the function

def runner(args):
    #Open the alignment file - expects BAM at the moment
    aln_file = pysam.AlignmentFile(args.input_file,'rb')
    
    #Make a dictionary of length = to the ref sequence
    coverage_dict = {}
    seq_length = aln_file.get_reference_length(args.ref_name)
    for i in range(0 ,seq_length):
        coverage_dict[i] = 0
    
    #Populate that dictionary with the read depth at each position
    for pileup_column in aln_file.pileup(args.ref_name, truncate=False, min_base_quality=0):
        coverage_dict[pileup_column.pos] = pileup_column.n
    
    #Create a list of lists containing the runs of positions below the "depth" value
    mask_pos_list = []
    position_list = []
    #The below could definitely be written more succinctly but it works...
    for pos in coverage_dict:
        if not args.inverse: # Added to give sites with at least X depth in bed format for samtools view
            if pos + 1 != seq_length: # Added to make sure last block to be masked is added to the file
                if coverage_dict[pos] < int(args.depth):
                    position_list.append(pos)
                else:
                    if position_list:
                        mask_pos_list.append(position_list)
                    position_list = []
            else:
                position_list.append(pos)
                mask_pos_list.append(position_list)
        else:
            if coverage_dict[pos] >= int(args.depth):
                position_list.append(pos)
            else:
                if position_list:
                    mask_pos_list.append(position_list)
                position_list = []
    
    #If there are no poisitons above the threshold
    if not mask_pos_list:
        no_cov = True
        if not args.inverse:
            mask_pos_list.append([0,seq_length - 1]) # -1 to preserve function with the 0vs1 base bit below

    #Write a "bcftools consensus" friendly file of regions to mask
    if mask_pos_list:
        mask_file = open(args.output_name + '.tsv', 'w')
        for mask_region in mask_pos_list:
            #bcftools expects your mask file to be one based so need to add one to all co-ordinates
            mask_file.write("%s\t%d\t%d\n" % (args.ref_name, mask_region[0] + 1, mask_region[-1] + 1 ))
        mask_file.close()

    #This bit actually masks a consensus file if you provide one 
    if args.fasta_to_mask:
        with open(args.fasta_to_mask) as cns: #Open the input
            for seq in SeqIO.parse(cns, "fasta"): #Get the seuqence bit
                seq_list = []
                if no_cov:
                    masked_seq = "N" * seq_length
                else:
                    for base in seq: #Next section extracts the sequence to a list to allow masking with the coordinates produced earlier
                        seq_list.append(base)
                    for line in mask_pos_list:
                        for i in line:
                            seq_list[i] = "N"
                    masked_seq = ''.join(seq_list)
                new_record = SeqRecord(Seq(masked_seq), id=seq.id, name=seq.name, description=seq.description)

                with open("%s.masked" % args.fasta_to_mask, "w") as output:
                    SeqIO.write(new_record, output, "fasta")



    aln_file.close()
    



def main():
    parser = argparse.ArgumentParser()
    '''
    Parse the command line arguments.
    '''
    parser = argparse.ArgumentParser(description='Creates a coverage mask to apply to your lovely consensus fasta.')

    
    optional_group = parser.add_argument_group('Optional')
    optional_group.add_argument('-d', '--depth', dest='depth', default="20",
                            help='If coverage is below this it will be masked. Default = 20')
    optional_group.add_argument('-r', '--ref-name', dest='ref_name', default="MN908947.3",
                            help='Name of ref the bam files were aligned to. Default = "MN908947.3"')
    optional_group.add_argument('-o', '--output-name', dest='output_name', default="depth_mask",
                            help='Prefix for the output. Default = "depth_mask"')
    optional_group.add_argument('-m', '--mask', dest='fasta_to_mask',
                            help='Mask a consensus sequence with your newly produced mask')
    optional_group.add_argument('-i', '--inverse', dest='inverse', action='store_true',
                            help='Return bed file of positions EQUAL OR ABOVE the chosen depth')
    optional_group.add_argument('-v', '--version', action='version', version='maskara 1.1.3',
                                help="Return Maskara version")

    parser.add_argument('input_file',
                            help='Path to the BAM file you want to create a mask for')


    args = parser.parse_args()
    runner(args)    
    
    
    
      
if __name__ == "__main__":
    main()

