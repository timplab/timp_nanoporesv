##From oxford website:

def barcode_extract(tary, prefix):
    ##HDF format
    import h5py
    ##Os IO stuff
    import os
    ##Tarball 
    import tarfile
    ##File recognize
    import glob
    ##shell util
    import shutil
    ##Random number
    import random
    ##Gzip lib
    import gzip

    
    ##Random number dir location for tar
    tmpdir='/tmp/'+str(random.randint(1e6, 9e6))
    
    os.mkdir(tmpdir)

    tarball=tarfile.open(name=tary, mode='r')

    tarball.extractall(path=tmpdir)

    filelist=glob.glob(tmpdir+'/*fast5')

    for filename in filelist:
        barcodesummary = '/Analyses/Barcoding_000/Summary/barcoding'
        barcodefastq = '/Analyses/Barcoding_000/Barcoding/Fastq'

        ##Load file
        hdf = h5py.File(filename, 'r')

        barcodeid = hdf[barcodesummary].attrs.get('barcode_arrangement')
        barcodefastq = hdf[barcodefastq]
                
        fqout=gzip.open(prefix+'_'+barcodeid+'.fastq.gz', 'ab')

        fqout.write(barcodefastq[()]+ '\n')

        fqout.close()
            
        hdf.close()        

    shutil.rmtree(tmpdir)
    

#barcode_extract('040615_ANbarcode_pass.tar.bz2', '040615_ANbarcode')
