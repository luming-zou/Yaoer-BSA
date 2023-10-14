#!/usr/bin/env python3

'''
MIT License

Copyright (c) 2017 Grok

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
'''

def process_vcf_file(vcf_filename):
    '''
    Process a VCF file and extract specific information.

    Parameters:
    - vcf_filename (str): The name of the VCF file to process.

    Returns:
    - None
    '''

    # Open the VCF file
    with open(vcf_filename) as vcf:
        # Iterate through each line in the VCF file
        for line in vcf:
            # Skip lines starting with '#', which are comments or header
            if not line.startswith('#'): 
                # Split the line into fields using tab as the delimiter
                CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT, S01, S02, S03, S04 = line.split('\t')
                
                # Extract information from the FORMAT field for each sample
                S01_d = S01.split(':')[1].split(',')
                S02_d = S02.split(':')[1].split(',')
                S03_d = S03.split(':')[1].split(',')
                S04_d = S04.split(':')[1].split(',')
                
                # Print the extracted information in a formatted way
                print('\t'.join(([CHROM, POS, REF, ALT] + S01_d + S02_d + S03_d + S04_d)))

# Example usage:
process_vcf_file("SNP_biallele.recode.vcf")
