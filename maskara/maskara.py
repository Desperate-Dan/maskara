#!/usr/bin/env python3

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import argparse
import pysam
import sys

#Main body of the function
def ref_finder(aln_file,reads):
    ref_stats = aln_file.get_index_statistics()
    refs_above_threshold = {}
    for ref in ref_stats:
        if ref[1] > int(reads):
            refs_above_threshold[ref[0]] = ref[1]
    return refs_above_threshold

def runner(args):
    #Open the alignment file - expects BAM at the moment
    aln_file = pysam.AlignmentFile(args.input_file,'rb')

    if args.mmm:
        refs_above_threshold = ref_finder(aln_file,args.reads)
        if not refs_above_threshold:
            print (f"No refs with more than {args.reads} reads!")
            sys.exit(0)
    else:
        refs_above_threshold = [args.ref_name]
    
    for ref in refs_above_threshold:
        #Make a dictionary of length = to the ref sequence
        coverage_dict = {}
        seq_length = aln_file.get_reference_length(ref)
        for i in range(0 ,seq_length):
            coverage_dict[i] = 0
        
        #Populate that dictionary with the read depth at each position
        for pileup_column in aln_file.pileup(ref, truncate=False, min_base_quality=0):
            #It turns out that using pileup_column.n ignores the min_base_quality arg so either need to enact the filter differently or change the couting method.
            #Could use the logic from piranha as that does individual base counts which might be nice to have, but keeping it simple here for now.
            good_bases = 0
            qualities_list = pileup_column.get_query_qualities()
            for base_q in qualities_list:
                if base_q > int(args.quality):
                    good_bases += 1
            coverage_dict[pileup_column.pos] = good_bases
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
            mask_file = open(ref + '_mask.tsv', 'w')
            for mask_region in mask_pos_list:
                #bcftools expects your mask file to be one based so need to add one to all co-ordinates
                mask_file.write("%s\t%d\t%d\n" % (ref, mask_region[0] + 1, mask_region[-1] + 1 ))
            mask_file.close()

    #This bit actually masks a consensus file if you provide one 
    if args.fasta_to_mask:
        if not args.mmm:
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
        else:
            print (f"Can't mask fastas with Multi-Map Mode, rerun for the respective refs individually. The refs with greater than {args.reads} reads mapped are:\n{refs_above_threshold}")


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
    optional_group.add_argument('-q', '--quality', dest='quality', default="20",
                            help='Choose the minimum base quality for consideration in coverage counting. Default = 20')
    optional_group.add_argument('--mmm', dest='mmm', action='store_true',
                                help="Multi-Map Mode: Get a depth mask for all references in your bam file with at least X reads")
    optional_group.add_argument('--reads', dest='reads', default="50",
                                help="Set the read limit for Multi-Map Mode. If a reference has at least this many reads it will have a depth mask made (note this is not the same as depth)")
    optional_group.add_argument('-v', '--version', action='version', version='maskara 1.1.7',
                                help="Return Maskara version")
    

    parser.add_argument('input_file',
                            help='Path to the BAM file you want to create a mask for')


    args = parser.parse_args()
    runner(args)    
    
    
    
      
if __name__ == "__main__":
    main()

