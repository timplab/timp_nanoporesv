# This is a python script for determining the number of reads which overlap a breakpoint region, as defined by a bed file
def split_bam_bed_overlap(bamfile):
    ##Support for gz files
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

    bammy=pysam.Samfile(bamfile, "rb")
    
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
    
    
    numreg=len(bedreg.index)

    comby=[]

    myfile=open(bamfile.replace('bam', 'cov'), 'wb')
    
    f = open(bamfile.replace('bam', 'pair.txt'), 'wb') 
    writer = csv.writer(f, delimiter = '\t')

        
    curread=""
    nonaligned=0
    aligned=0
    numsplits=[0]*6

    for read in bammy.fetch(until_eof=True):    
        if curread != read.qname: 
            curread=read.qname
            if read.reference_id == -1:
                nonaligned += 1
            else:
                aligned += 1                             
                comby.append(str(bedreg.name[bedreg.coinc.nonzero()[0]].tolist()))                
                bedreg.coinc = 0
    
        if read.reference_id != -1:
            
            #print (read.reference_end-read.reference_start)
            
            for idx, reg in bedreg.iterrows():
                
                ##chr must match AND both coordinates must not be on same side
                if ( (bammy.getrname(read.reference_id)==reg['chr']) and 
                     (not(((read.reference_start < reg['start']) and (read.reference_end < reg['start'])) or 
                          ((read.reference_start > reg['end']) and (read.reference_end > reg['end']))))):                    
                    ##If coincident and running tally of how many per region
                    bedreg.loc[idx,'coinc']=1
                    bedreg.loc[idx,'overlaps'] += 1

                    ##How far off upstream and downstream from target was alignment on avg and variance
                    ## Knuth's running variance/mean algo
                    deltaup=((reg['start']-read.reference_start)-bedreg.loc[idx,'upstream_mean'])
                    deltadown=((read.reference_end-reg['end'])-bedreg.loc[idx,'downstream_mean'])
                    bedreg.loc[idx,'upstream_mean'] += deltaup/bedreg.loc[idx,'overlaps']
                    bedreg.loc[idx,'downstream_mean'] += deltadown/bedreg.loc[idx,'overlaps']
                    
                    bedreg.loc[idx,'upstream_var'] += deltaup*((reg['start']-read.reference_start)-bedreg.loc[idx,'upstream_mean'])
                    bedreg.loc[idx,'downstream_var'] += deltadown*((read.reference_end-reg['end'])-bedreg.loc[idx,'downstream_mean'])

                    ##How many perfect alignments
                    if (reg['start']==read.reference_start):
                        bedreg.loc[idx,'upperfect'] += 1
                    if (reg['end']==read.reference_end):
                        bedreg.loc[idx,'downperfect'] += 1

                        
    
    ##on last one do it too
    comby.append(str(bedreg.name[bedreg.coinc.nonzero()[0]].tolist()))
    bedreg.coinc = 0

    bammy.close()

    myfile.write('Aligned: ' + str(aligned) + ' Nonaligned: ' + str(nonaligned) + ' Total: ' + str(nonaligned+aligned)+'\n')
    
    #print bedreg

    ##Finish variance algo
    bedreg['upstream_var'] /= bedreg['overlaps']-1
    bedreg['downstream_var'] /= bedreg['overlaps']-1

    bedreg['upstream_std']=np.sqrt(bedreg['upstream_var'])
    bedreg['downstream_std']=np.sqrt(bedreg['downstream_var'])
    
    freqs=Counter(comby)
    #print freqs

    for key, count in freqs.most_common():
        writer.writerow([key] + [count])

    #print comby

    bedreg.to_csv(myfile, sep="\t")
        
    myfile.close()
    f.close()

import sys

split_bam_bed_overlap(sys.argv[1])



