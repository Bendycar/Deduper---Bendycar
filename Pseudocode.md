DEFINING THE PROBLEM: 
PCR is a necessary part of the RNA-seq workflow in order to get a sufficient signal to make sure we sequence all expressed RNA. However, since we are interested in the abundance of our RNA data, we need to eliminate these PCR duplicates. Because PCR creates replicas unevenly across our samples, this is not just as simple as dividing by the amount duplicated. This problem is also complicated by the fact that we may have legitimate biological duplicate reads that we do not want to remove, again because we are interested in abundance. To solve this, Unique Molecular Identifiers are used to differentiate between duplicates of the same exact molecule, and duplicates of the same position on the genome.

A read is considered to be a PCR duplicate if it shares ALL of the following with another read:

1. Chromosome 
2. Position (Adjusted for soft clipping)
3. Strandedness
4. UMI

In this case, we want to save only the first instance of a PCR duplicate to our output file.

To check if any given read shares these characteristics with another, I will create a set of [Position, Strandedness, and UMI]. To save memory, I will clear this set between every chromosome (the files will be sorted by SAMTools before running). For this same reason, I will not be storing the chromosome data, as we can assume every item in the set has the same chromosome. Because we still want to write out the first instance of a PCR duplicate, I will write out any read that is not already in this set, except for those with invalid UMIs. 

---------------------------------------------------------------------------------------------------------------------------------------------------
FUNCTIONS:

get_position(line, strandedness):
```Takes in a single SAM file line (already split by tab), returns the 5' starting position, ADJUSTED FOR SOFT CLIPPING```
#Must be able to handle minus strands as well -- ADD all D's, M's, and S's to the RIGHT of M to position. 
cigar_string = pull out whichever position this is in
position = wherever that is in split SAM file
soft_clip_adjustment = the number proceeding S if it is before the M, otherwise returns 0
return position - soft_clip_adjustment

Example Input: get_position('NS500451:154:HWKTMBGXX:1:11101:25071:87858:TAGCAAGG	0	2	76718924	36	71M	*	0	0	GTTCCTCTTGCATTAATGGCTGTGACCCGATGTTCATATTCTAAGCCTTCAATAAGGCCGGTGGAGCGGTA	6/EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEAEEEE	MD:Z:71	NH:i:1	HI:i:1	NM:i:0	SM:i:36	XQ:i:40	X2:i:0	XO:Z:UU
')
Example Output: 76718924

get_chr(line):
```Takes in a single SAM file line (already split by tab), returns the chromosome number```
return chromosome

Example Input: get_chr('NS500451:154:HWKTMBGXX:1:11101:25071:87858:TAGCAAGG	0	2	76718924	36	71M	*	0	0	GTTCCTCTTGCATTAATGGCTGTGACCCGATGTTCATATTCTAAGCCTTCAATAAGGCCGGTGGAGCGGTA	6/EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEAEEEE	MD:Z:71	NH:i:1	HI:i:1	NM:i:0	SM:i:36	XQ:i:40	X2:i:0	XO:Z:UU
')
Example Output: 2

get_strandedness(line):
```Takes in a single SAM file line (already split by tab), returns + or - for strand```
Parse bitwise flag to determine strandedness
Return strandedness (+ or -)

Example Input: get_strandedness(''NS500451:154:HWKTMBGXX:1:11101:25071:87858:TAGCAAGG	0	2	76718924	36	71M	*	0	0	GTTCCTCTTGCATTAATGGCTGTGACCCGATGTTCATATTCTAAGCCTTCAATAAGGCCGGTGGAGCGGTA	6/EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEAEEEE	MD:Z:71	NH:i:1	HI:i:1	NM:i:0	SM:i:36	XQ:i:40	X2:i:0	XO:Z:UU
')
Example Output: +


get_UMI(line):
```Takes in a single SAM file line (already split by tab), extracts the UMI from header```
Split header to extract UMI
Return UMI

Example Input: get_UMI('NS500451:154:HWKTMBGXX:1:11101:25071:87858:TAGCAAGG	0	2	76718924	36	71M	*	0	0	GTTCCTCTTGCATTAATGGCTGTGACCCGATGTTCATATTCTAAGCCTTCAATAAGGCCGGTGGAGCGGTA	6/EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEAEEEE	MD:Z:71	NH:i:1	HI:i:1	NM:i:0	SM:i:36	XQ:i:40	X2:i:0	XO:Z:UU
')
Example Output: TAGCAAGG

def get_args(): 
    parser.add_argument("-f", "--input", help="Absolute file path to sorted SAM file.", type=str, required=True)
    parser.add_argument("-o", "--outfile", help="Absolute file path to deduplicated SAM file.", type=str, required=True)
    parser.add_argument("-u", "--umi", help="File containing the list of UMIs.", type=str, required=True)
    parser.add_argument("-h", "--help", help="Useful information about the program.", type=str, required=False)
    return parser.parse_args()

Code should run just by: ./Carr_deduper.py -u STL96.txt -f <in.sam> -o <out.sam>

---------------------------------------------------------------------------------------------------------------------------------------------------

PSEUDOCODE:

open SAM_IN and SAM_OUT
UMI_set = set(UMIs)
RC_UMI_set = rev_comp(UMI_set) #Need to account for the fact the the minus strands have reverse complemented the UMI. This should be the only case where it matters, because all UMIs of reads on the minus strand will be compared to the RC UMI
seen = {}
chrom_num = 1

If needed write out the first few header lines until we hit the actual reads (maybe check while line starts with @, but could be risky...)
```
while True:
    line = readline(SAM_IN)
    if line == "":
        break

    spline = line.split['\t']
    chr = get_chr(spline)
    position = get_position(spline)
    strandedness = get_strandedness(line)
    UMI = get_UMI(line)

    if UMI not in UMI_set or RC_UMI_set if minus strand:
        continue
    
    if chr != chrom_num: 
        clear 'seen' set #This ensures that we are only storing one chromosome worth of data in memory at a time
        chrom_num = chr
    
    identifiers = [position, strandedness, UMI]

    if identifiers not in seen: #On the first read of a new chromosome, we will have just cleared the set so this will alwyas be true
        add to seen
        outfile.write(line)
```
---------------------------------------------------------------------------------------------------------------------------------------------------
TESTCASES:

I believe we will need to test the following cases: 

1. Unique read. Line 1 in test file, every element different from all others DONE
2. Trivial duplicates. Lines 2 and 3 in my test file, just the same read copied and pasted DONE
3. Adjusted position duplicates. Lines 4 and 5 in my test file. Everything same except for position, which will become the same once adjusted for soft clipping. DONE

2 and 3 should be the only cases where duplicates are removed -- the remaining cases are "almost duplicates" that we want to keep.

4. Biological duplicates (everything but UMI). Lines 6 and 7 in test file. Exact same information EXCEPT UMI DONE
5. Unadjusted position duplicates. Lines 8 and 9. Same position without soft clipping adjustment, after which position is different DONE
6. Everything but position duplicates. Lines 10 and 11. Everything is same except unadjusted position DONE
7. Everything but adjusted position duplicates. Lines 12 and 13. Everything is same except position after soft clipping adjustment (not necessarily same initially, like in case 4). DONE
8. Everything but chromosome duplicates. Lines 14 and 15. DONE
9. Everything but strandedness duplicates. Lines 16 and 17. DONE
10. Incorrect adjusted position duplicates. Lines 18 and 19. Two reads that would be adjusted position duplicates if the program incorrectly adjusts position when the soft clipping is to the RIGHT of the 'M' DONE
11. Invalid UMI, line 20 DONE
12. Invalid UMI on minus strand, line 21 DONE

Which should correspond to the following lines in my output file:
1. Unique read DONE
2. One instance of the trivial duplicates DONE
3. One instance of adjusted position duplicates DONE
4. First instance of biological duplicates DONE
5. Second instance of biological duplicates DONE
6. First instance of unadjusted position duplicates DONE
7. Second instance of unadjusted position duplicates DONE
8. First instance of everything but position duplicates DONE
9. Second instance of everything but position duplicates DONE
10. First instance of everything but adjusted position duplicates DONE
11. Second instance of everything but adjusted position duplicates DONE
12. First instance of everything but chromosome duplicates DONE
13. Second instance of everything but chromosome duplicates DONE
14. First instance of everything but strandedness duplicates DONE
15. Second instance of everything but strandedness duplicates DONE
16. First instance of incorrect adjusted position duplicates DONE
17. Second instance of incorrect adjusted position duplicates DONE
