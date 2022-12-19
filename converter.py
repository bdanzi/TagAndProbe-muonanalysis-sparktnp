import os
import subprocess
import glob
import getpass
from pyspark.sql import SparkSession

from registry import registry
from dataset_allowed_definitions import get_allowed_sub_eras


def run_convert(spark, particle, resonance, era, dataTier, subEra, customDir='', baseDir='', use_pog_space=False, use_local=False):
    '''
    Converts a directory of root files into parquet
    '''

    if baseDir == '':
        if use_pog_space is False:
            inDir = os.path.join('hdfs://analytix/user', getpass.getuser())
            outDir = os.path.join('hdfs://analytix/user', getpass.getuser())
        else:
            inDir = 'hdfs://analytix/cms/muon_pog'
            outDir = 'hdfs://analytix/cms/muon_pog'
        inDir = os.path.join(inDir, 'root', customDir, particle, resonance, era, dataTier, subEra)
        outDir = os.path.join(outDir, 'parquet', customDir, particle, resonance, era, dataTier, subEra)
    else:
        if use_local is False and 'hdfs' not in baseDir:
            print('>>>>>>>>> Warning! Custom baseDir given to convert but `useLocalSpark` flag not enabled and no `hdfs` in baseDir.')
            print('>>>>>>>>> Distributed spark clusters can only read files in hdfs. Make sure baseDir is an `hdfs` path, e.g.:')
            print('>>>>>>>>> `hdfs://analytix/user/[user]/[your-custom-dir]`')
            print('>>>>>>>>> Or else, if this is a local test, use the `--useLocalSpark` flag (see command help).')
        inDir = baseDir
        outDir = baseDir

    # The following glob command works only on the edge nodes, which has a fuse-style hdfs mountpoint    
    # if 'hdfs' in inDir:
    #     inDir = inDir.replace('hdfs://analytix', '/hdfs/analytix.cern.ch')
    #     fnames = glob.glob(os.path.join(inDir, f'{inDir}/*.root'))
    # else:
    #     fnames = glob.glob(f'{inDir}/*.root')
    # fnames = [f.replace('/hdfs/analytix.cern.ch', 'hdfs://analytix') for f in fnames]
    
    # Make sure path is in hdfs format (not fuse-style format)
    inDir = inDir.replace('/hdfs/analytix.cern.ch', 'hdfs://analytix')
    outDir = outDir.replace('/hdfs/analytix.cern.ch', 'hdfs://analytix')

    # The following glob command works in both lxplus and edge nodes
    cmd = "hdfs dfs -find {} -name '*.root'".format(inDir)
    fnames = subprocess.check_output(cmd, shell=True).strip().split(b'\n')
    fnames = [fname.decode('ascii') for fname in fnames]
    
    outName = os.path.join(outDir, 'tnp.parquet')
    if use_local is True and 'hdfs' not in outName:
        outName = 'file://' + outName
    else:
        outName = outName.replace('/hdfs/analytix.cern.ch', 'hdfs://analytix') # just in case

    print(f'>>>>>>>>> Path to input root files: {inDir}')
    print(f'>>>>>>>>> Path to output parquet files: {outName}')
    
    # treename = 'tpTree/fitter_tree' # old tnp tool tree name
    # treename = 'muon/tree' # old miniAOD tree name
    treename = 'muon/StandAloneEvents'

    print(f'>>>>>>>>> Number of files to process: {len(fnames)}')
    if len(fnames) == 0:
        print('>>>>>>>>> Error! No ROOT files found to convert with desired options.')
        print('>>>>>>>>> Exiting...')
        return
    print(f'>>>>>>>>> First file: {fnames[0]}')

    # process batchsize files at a time
    batchsize = 100
    new = True
    while fnames:
        current = fnames[:batchsize]
        fnames = fnames[batchsize:]

        rootfiles = spark.read.format("root")\
                         .option('tree', treename)\
                         .load(current)
                         
        # merge rootfiles. chosen to make files of 8-32 MB (input)
        # become at most 1 GB (parquet recommendation)
        # https://parquet.apache.org/documentation/latest/
        # .coalesce(int(len(current)/32)) \
        # but it is too slow for now, maybe try again later
        if new:
            rootfiles.write.parquet(outName)
            new = False
        else:
            rootfiles.write.mode('append')\
                     .parquet(outName)


def run_all(particle, resonance, era, dataTier, subEra=None, customDir='', baseDir='', use_pog_space=False, use_local=False):

    if subEra is not None:
        subEras = [subEra]
    else:
        subEras = get_allowed_sub_eras(resonance, era)
        # subEras by default includes whole era too, so remove for convert
        subEras.remove(era)

    local_jars = ','.join([
        './laurelin-1.0.0.jar',
        './log4j-api-2.13.0.jar',
        './log4j-core-2.13.0.jar',
    ])
    
    spark = SparkSession\
        .builder\
        .appName("TnP")\
        .config("spark.jars", local_jars)\
        .config("spark.driver.extraClassPath", local_jars)\
        .config("spark.executor.extraClassPath", local_jars)\
        .config("spark.dynamicAllocation.maxExecutors", "100")\
        .config("spark.driver.memory", "10g")\
        .config("spark.executor.memory", "10g")\
        .config("spark.sql.shuffle.partitions", "500")\
        .config("spark.executor.cores", "1")

    
    if use_local is True:
        spark = spark.master("local")

    spark = spark.getOrCreate()

    print('\n\n------------------ DEBUG ----------------')
    sc = spark.sparkContext
    print(sc.getConf().toDebugString())
    print('---------------- END DEBUG ----------------\n\n')

    for subEra in subEras:
        print('\n>>>>>>>>> Converting:', particle, resonance, era, subEra)
        run_convert(spark, particle, resonance, era, dataTier, subEra, customDir, baseDir, use_pog_space, use_local)

    spark.stop()
