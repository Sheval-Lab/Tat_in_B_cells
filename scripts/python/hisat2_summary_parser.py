#!/usr/bin/env python
# coding: utf-8

# In[1]:


import os


# In[19]:


# Input directory with summary files
file_path = '../../data/summaries/'

# Output file
out_file = 'hisat2_summary.tsv'
out = open(str(file_path+out_file), 'w')
# Write header into output file
out.write('\t'.join(['sample', 'total_reads', 'zero_al', 'one_al', 'mult_al\n']))


# List files in the input dir
files = os.listdir(file_path)


for f in files:
    # Extract sample name from file name
    sample = f.strip('.log')
    # Read each summary file
    file = open(file_path+f, 'r')
    # Process file only if it contains 'log' (format suffix)
    if 'log' in f:
#         print(str(file_path+f))
        k = 0
        for line in file:
            k += 1
            # Line with total reads number
            if k == 1:
                total_reads = int(line.split()[0])
            # Line with number of reads failed to align
            if k == 3:
                zero_al = int(line.split()[0])
            # Line with uniquely aligned reads number
            if k == 4:
                one_al = int(line.split()[0])
            # Line with multi-aligned reads number
            if k ==5:
                mult_al = int(line.split()[0])
        # Write a line into tsv file for each summary file
        out.write('\t'.join([sample, str(total_reads), str(zero_al), str(one_al), str(mult_al)+'\n']))
    file.close()
        
out.close()

