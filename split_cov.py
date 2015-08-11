# This is a python script for determining the number of reads which overlap a breakpoint region, as defined by a bed file
def split_bam_bed_overlap(bamfile):
    
    ##Import necessary modules
    import pysam
    import string
    import gzip
    import csv
    from Bio import SeqIO

    import pysam
    import itertools

    import numpy as np
    import pandas as pd

    from collections import Counter

    ##Load the bam file
    bammy=pysam.Samfile(bamfile, "rb")
    
    ##Load the list of expected regions and initialize columns for statistical information
    bedreg = pd.read_csv('p16_amp.bed', sep='\t', names=('chr', 'start', 'end', 'name'))
    bedreg['coinc']=0
    bedreg['overlaps']=0
    bedreg['upperfect']=0
    bedreg['upstream_mean']=0.0
    bedreg['upstream_std']=0.0
    bedreg['downperfect']=0
    bedreg['downstream_mean']=0.0
    bedreg['downstream_std']=0.0
    bedreg['upstream_var']=0.0
    bedreg['downstream_var']=0.0
    
    
    ##Open files and initailize variables to be used subsequently in for loops
    numreg=len(bedreg.index)

    comby=[]

    myfile=open(bamfile.replace('bam', 'cov'), 'wb')
    
    f = open(bamfile.replace('bam', 'pair.txt'), 'wb') 
    writer = csv.writer(f, delimiter = '\t')

        
    curread=""
    nonaligned=0
    aligned=0
    numsplits=[0]*6

    ##Loop through each read in the bam file
    for read in bammy.fetch(until_eof=True):    
        ##Each read may align multiple times, resulting in multiple reads with the same qname in the bam file. 
        if curread != read.qname: ##If the read changes, record whether it aligned or not.
            curread=read.qname
            if read.reference_id == -1:
                nonaligned += 1
            else:
                aligned += 1 
                ##Keep track of which regions are overlapped by each read. bedreg.coinc is set to 1 for each region a read touches
                comby.append(str(bedreg.name[bedreg.coinc.nonzero()[0]].tolist()))                
                bedreg.coinc = 0
        
        ##Check every read that aligned against every region
        if read.reference_id != -1:
                        
            for idx, reg in bedreg.iterrows():
                
                ##Check that the read overlaps the region. chr must match AND both coordinates must not be on same side
                if ( (bammy.getrname(read.reference_id)==reg['chr']) and 
                     (not(((read.reference_start < reg['start']) and (read.reference_end < reg['start'])) or 
                          ((read.reference_start > reg['end']) and (read.reference_end > reg['end']))))):                    
                    
                    ##Keep a running tally of how many reads aligned per region.
                    bedreg.loc[idx,'coinc']=1
                    bedreg.loc[idx,'overlaps'] += 1

                    ##Keep track of how far off upstream and downstream from target was alignment on avg and variance
                    ##Use Knuth's running variance/mean algo
                    deltaup=((reg['start']-read.reference_start)-bedreg.loc[idx,'upstream_mean'])
                    deltadown=((read.reference_end-reg['end'])-bedreg.loc[idx,'downstream_mean'])
                    bedreg.loc[idx,'upstream_mean'] += deltaup/bedreg.loc[idx,'overlaps']
                    bedreg.loc[idx,'downstream_mean'] += deltadown/bedreg.loc[idx,'overlaps']
                    
                    bedreg.loc[idx,'upstream_var'] += deltaup*((reg['start']-read.reference_start)-bedreg.loc[idx,'upstream_mean'])
                    bedreg.loc[idx,'downstream_var'] += deltadown*((read.reference_end-reg['end'])-bedreg.loc[idx,'downstream_mean'])

                    ##Keep a running tally of how may reads aligned perfectly
                    if (reg['start']==read.reference_start):
                        bedreg.loc[idx,'upperfect'] += 1
                    if (reg['end']==read.reference_end):
                        bedreg.loc[idx,'downperfect'] += 1

                        
    
    ##Iterate one last time since fetch attribute stops at the last read. 
    comby.append(str(bedreg.name[bedreg.coinc.nonzero()[0]].tolist()))
    bedreg.coinc = 0

    bammy.close()
    
    ##Print alignment information as the first line of the output csv file.
    myfile.write('Aligned: ' + str(aligned) + ' Nonaligned: ' + str(nonaligned) + ' Total: ' + str(nonaligned+aligned)+'\n')
    
    
    ##Finish variance algo
    bedreg['upstream_var'] /= bedreg['overlaps']-1
    bedreg['downstream_var'] /= bedreg['overlaps']-1

    bedreg['upstream_std']=np.sqrt(bedreg['upstream_var'])
    bedreg['downstream_std']=np.sqrt(bedreg['downstream_var'])
    
    
    freqs=Counter(comby)
    #print freqs
    
    ##Print which combinations of regions is most common.
    for key, count in freqs.most_common():
        writer.writerow([key] + [count])

    ##Write the region, alignment, and up/downstream information to a csv. 
    bedreg.to_csv(myfile, sep="\t")
        
    myfile.close()
    f.close()

import sys

split_bam_bed_overlap(sys.argv[1])



