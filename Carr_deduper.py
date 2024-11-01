#!/usr/bin/env python

import argparse
import re

def get_args(): 
    parser = argparse.ArgumentParser(description="Python program for reference based PCR duplicate removal. Input SAM file must be presorted with SAMTools, and a list of valid UMIs must be supplied.")
    parser.add_argument("-f", "--input", help="Absolute file path to sorted SAM file.", type=str, required=True)
    parser.add_argument("-o", "--outfile", help="Absolute file path to deduplicated SAM file.", type=str, required=True)
    parser.add_argument("-u", "--umi", help="File containing the list of UMIs.", type=str, required=True)
    return parser.parse_args()

args = get_args()
input_file = args.input
output_file = args.outfile
UMI_file = args.umi

def get_umi_set(UMI_file: str) -> set:
    '''Takes in the list of UMI sequences provided by the -u flag, converts to a set for checking valid UMIs in every read. Input file must contain exactly one valid UMI sequence per line'''
    UMI_set = set()
    with open(UMI_file, 'r') as UMI_in:
        while True:
            line = UMI_in.readline().strip()
            if line == '':
                break
            UMI_set.add(line)
    return UMI_set

def get_line_info(line: str) -> tuple:
    '''Takes in unmodified (unsplit) SAM file line, returns Chromosome, UMI, Strandedness (+ or -), and 5' starting position adjusted for soft clipping'''
    spline = line.split()
    chr = spline[2] # Never casted to int so we can handle potential X, Y, or M chromosomes
    umi = spline[0].split(':')[-1] # Example header line: '18:NS500451:154:HWKTMBGXX:1:11101:69992:67325:CTAGAGGA'. The UMI is found at the very end.
    cigar = spline[5]
    pos = int(spline[3])
    strandedness = get_strandedness(int(spline[1])) # Uses the bitwise flag in the SAM file at position 1

    pos = adjust_position(pos, strandedness,cigar) # Adjustment for softclipping and strandedness

    return (chr, umi, strandedness, pos)
    
def get_strandedness(bitwise_flag: int) -> str:
    '''Takes in bitwise flag from SAM file, uses bitwise comparison to determine the strandedness of the read'''
    if bitwise_flag & 16 == 16: # Bitwise flag 16s place is defined as "SEQ being reverse complemented"
        return '-' # A reverse complemented string implies that it is considered the minus strand
    else:
        return '+' 

def adjust_position(position: int, strandedness: str, cigar: str) -> int:
    '''The SAM file contains leftmost 1-based starting position, which may differ from the 5' starting position, which we need to use to determine duplication.
    In a read on the plus strand, the only adjustment that needs to be made is correcting for soft clipping on the 5' end. To correct for this, that value is subtracted from the leftmost starting position.
    In the case of a read on the minus strand, every N, M, and D in the CIGAR string must be added to the leftmost starting position value, as well as soft clipping on the 5' end.
    To get the true 5' starting position, technically we should subtract 1 from the position to account for being 1-based, but this is unnecessarily as long as we are consistent with every read.'''
    if strandedness == '+':
        softclip = re.match(r'([0-9]+)S', cigar) # Regex match only returns patterns at the START of a string that match the expression
        if softclip: # Therefore if there is no S, or the S is at the end of the cigar string, this expression will return None
            softclip_adjustment = softclip.group(1) # Returns the first captured group, which will be the number in front of the leading S
            position -= int(softclip_adjustment)
    else:
        consumes_reference = re.findall(r'([0-9]+)[NDM]', cigar) # Findall will return a list of strings matching the captured portion of my pattern -- in this case, just the numbers associated with elements that consume reference in cigar string
        for i in consumes_reference:
            position += int(i)
        softclip = re.search(r'([0-9]+)S$', cigar)# The $ ensures this pattern will only match at the END of the string, meaning 5' softclipping that we need to add to our position
        if softclip:
            softclip_adjustment = softclip.group(1)
            position += int(softclip_adjustment)
    
    return position

def main():
    with open(output_file, 'w') as outfile, open(input_file, 'r') as infile:
        UMI_set = get_umi_set(UMI_file)
        seen = set()
        invalid_UMIs = 0 # Variables for answering final questionnaire
        duplicates_removed = 0
        
        while True: # First loop writes out all header lines
            line = infile.readline()
            if not line.startswith('@'): #First actual read, need to store this one before moving on to next loop
                chr, umi, strandedness, pos = get_line_info(line)
                current_chr = chr
                if umi in UMI_set:
                    information = (umi,strandedness,pos) #Don't need to include chromosome because this set will be wiped every time we reach a new chromosome
                    seen.add(information)
                    outfile.write(line)
                else: # Invalid UMI in first read
                    invalid_UMIs +=1
                break
            else:
                outfile.write(line) # Writing the header lines

        while True: # Loop to actually deduplicate
            line = infile.readline()
            if line == '':
                break
            chr, umi, strandedness, pos = get_line_info(line)
            if umi not in UMI_set: # We never want to do anything with an invalid UMI
                invalid_UMIs += 1 
                continue
            if chr != current_chr: # This conditional wipes the set every time we come across a new chromosome for the sake of memory efficiency
                seen.clear()
                current_chr = chr
            information = (umi,strandedness,pos)
            if information not in seen:
                seen.add(information)
                outfile.write(line)
            else:
                duplicates_removed += 1
    
    print(f"Number of reads removed due to invalid UMI sequences: {invalid_UMIs}")
    print(f"Number of reads removed due to PCR duplication: {duplicates_removed}")

if __name__ == '__main__':
    main()