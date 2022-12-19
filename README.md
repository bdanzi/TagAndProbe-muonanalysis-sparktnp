# spark-tnp
Tag and probe analysis using Apache Spark.

## Before you begin
This package uses Apache Spark clusters.
More details on CERN's Apache Spark can be found [here](https://hadoop-user-guide.web.cern.ch/hadoop-user-guide/spark/Using_Spark_on_Hadoop.html).

**Important:** If you want to use the CERN analytix cluster (which is much faster to startup than the k8s cluster),
you need to request access to the cluster in the help document found [here](https://hadoop-user-guide.web.cern.ch/hadoop-user-guide/getstart/access.html).

## Quick start

The following will produce a set of example efficiencies (assuming you can run on analytix) (assuming Git version > 2.13):

```bash
git clone https://gitlab.cern.ch/cms-muonPOG/spark_tnp.git
cd spark_tnp
source env.sh
kinit
./tnp_fitter.py flatten muon generalTracks Z Run2018_UL configs/muon_example.json --baseDir ./example
./tnp_fitter.py fit muon generalTracks Z Run2018_UL configs/muon_example.json --baseDir ./example
./tnp_fitter.py prepare muon generalTracks Z Run2018_UL configs/muon_example.json --baseDir ./example
```

## Interactive notebook

There are example notebooks in the [notebooks](notebooks) directory demonstrating the use of these tools interactively.
A good starting point is to follow the instructions in [MuonTnP.ipynb](notebooks/MuonTnP.ipynb).
These notebooks use [https://swan.cern.ch](https://swan.cern.ch) to connect to the Apache Spark clusters at CERN.

## Command line setup

There are a couple of ways you can run. Either connect to the hadoop cluster edge nodes or directly on lxplus.

**Update (January 2021):** It is highly recommended to use the edge nodes because they create fuse mountpoints for both eos and hdfs filesystems. This means one can move files easily from eos (e.g. flat ROOT ntuples produced with the [MuonAnalyzer package](https://gitlab.cern.ch/cms-muonPOG/muonanalysis-muonanalyzer/)) to hdfs using standard `cp` and `mv` commands, and tab completion is enabled. To move files on lxplus, dedicated `hdfs dfs` commands are needed instead (still doable, just not very convenient). Note that if you do work on the edge nodes you will need to briefly get back to lxplus to submit fitting condor jobs (see below), since the edge nodes do not set up condor.

The flattening jobs are run on spark clusters and the data is read from an hdfs cluster. The fitting jobs are still run on condor since they depend on ROOT libraries. The default (and preferred) way to run flattening jobs is to use the `analytix` spark and hdfs cluster.

### Edge node

Connect to the hadoop edge node (from within the CERN network, e.g. connect first to lxplus and then to the edge node or use a VPN):

```bash
ssh it-hadoop-client
```

Setup the environment:

```bash
kinit
source env.sh
```

### LXPLUS

Connect to LXPLUS:

```bash
ssh lxplus.cern.ch
```

Setup the environment:

```bash
kinit
source env.sh
```

**Note**: Do not forget to make sure you have a valid kerberos token with `kinit`.

### Recommended

To produce the new common schema jsons recommended by xPOG in addition to the classic style jsons, an additional library ([correctionlib]((https://github.com/nsmith-/correctionlib))) is needed:

```bash
pip install --user correctionlib==2.0.0rc4
```

### Optional

Install `tqdm` packaged for a nice progress bar.

```bash
pip install --user tqdm
```

## Tag-and-probe steps

The tag-and-probe process is broken down into several parts:

1. Creation of the flat ROOT tag-and-probe trees (not shown here, use the [MuonAnalyzer package](https://gitlab.cern.ch/cms-muonPOG/muonanalysis-muonanalyzer/))
2. Conversion of the ROOT TTree into the parquet data format (convert)
3. Reduce the data into binned histograms with spark (flatten)
4. Fit the resulting histograms (fit)
5. Extraction of efficiencies and scale factors (prepare)

These steps are controlled with the [tnp_fitter.py](tnp_fitter.py) script.
For help with the script run:
```bash
./tnp_fitter.py -h
```

The most important argument to pass is the configuration file
that controls what kind of fits are produced.
See detailed documentation in the [configs](configs) directory.

New tag-and-probe datasets will need to be registered in the [data](data) directory.

### 2. Conversion to parquet

The conversion to parquet vastly speeds up the later steps.
We will use [laurelin](https://github.com/spark-root/laurelin) to
read the root files and then write them in the parquet data format.
There are two possible approaches: using `k8s` and using `analytix`.

Conversion with `k8s` currently only works if you are using [https://swan.cern.ch](https://swan.cern.ch).
Use the [RootToParquet](notebooks/RootToParquet.ipynb) notebook as a guide.
The output should be written to `analytix`.

Conversion with `analytix` requires you to first copy your root files
to `hdfs://analytix`. There is an issue with reading root files from `eos`
on `analytix` that needs to be understood.

**Update (January 2021):** Using `analyix` is the recommended way to convert root files to parquet. The tips above to connect to the edge node and move files using the fuse mount makes the process much smoother.

The following should be executed when you are connected to the edge node:

```bash
cp /eos/user/[u]/[user]/[path-to-files]/*.root /hdfs/analytix.cern.ch/user/[user]/[path-to-out-dir]
```

or

```bash
cp /eos/cms/store/[path-to-files]/*.root /hdfs/analytix.cern.ch/cms/muon_pog/[path-to-out-dir]
```

depending on where your files are located and your access level.

If on lxplus, use `hdfs` commands to copy the files (note the slightly different syntax):

```bash
hdfs dfs -cp root://eoscms.cern.ch//eos/cms/store/[path-to-files]/*.root hdfs://analytix/[path-to-out-dir]
```

<details><summary>Optimizing ROOT file sizes to work with spark (click to expand) </summary>

**Update (May 2021):** Our ntuples have grown considerably in size with the latest developments. Hadding the files now can choke spark so keeping the output files separate is probably the safest way to go, even if not the absolute optimal in terms of speed. The instructions below are kept for reference only.

Spark works best with a few large-ish files (1 GB) as opposed to lots of small files (few MBs), which are produced by CRAB. Before converting to parquet, it is recommended to check the CRAB output file sizes of your ntuples and if they're very small to hadd them into a few large root files. A script `hadd_by_size.sh` exists in scripts/ to automatically convert a folder of small ROOT files into equisized 1 GB files. The syntax is:

```bash
./scripts/hadd_by_size.sh [path-to-crab-output]
```

This will produce several files named `haddOut_$num_$UUID.root` in the local directory, each approximately 1 GB in size. These files can then be transferred to `hdfs` to be converted into parquet. Note this step is not required but it is recommended for speed up gains later on.

</details>

You also will need to download the `jar` files to add
to the spark executors:

```bash
bash setup.sh
```

Once copied, you can use:

```bash
./tnp_fitter.py convert [particle] [resonance] [era] [[subera]]
```

By default the converter looks for root files in

```
/hdfs/analytix.cern.ch/user/[user]/root/[particle]/[resonance]/[era]/[subera]/
```

though a custom directory can be specified. The output parquet files will be placed in

```
/hdfs/analytix.cern.ch/user/[user]/parquet/[particle]/[resonance]/[era]/[subera]/tnp.parquet
``` 

Any new datasets must be registered in `registry.py` following instructions [here](data/README.md).

### 3. Flatten histograms with spark

This step uses the converted parquet data format to efficiently aggregate
the efficiency data into binned histograms.

```bash
./tnp_fitter.py flatten -h
```

For example, to flatten all histograms for the Run2017 Legacy muon scale factors from Z:

```bash
./tnp_fitter.py flatten muon generalTracks Z Run2017_UL configs/muon_pog_official_run2_Z_2017.json
```

You can optionally filter the efficiencies and shifts you flatten with the `--numerator`,
`--denominator`, and `--shiftType` arguments. Thus, to only flatten the nominal histograms do:
```bash
./tnp_fitter.py flatten muon generalTracks Z Run2017_UL configs/muon_pog_official_run2_Z_2017.json --shiftType Nominal
```

**Note:** running this on lxplus will give the following warnings:

>WARN NativeCodeLoader: Unable to load native-hadoop library for your platform... using builtin-java classes where applicable  
>WARN DomainSocketFactory: The short-circuit local reads feature cannot be used because libhadoop cannot be loaded.  
>WARN Client: Neither spark.yarn.jars nor spark.yarn.archive is set, falling back to uploading libraries under SPARK_HOME.  
>WARN YarnSchedulerBackend$YarnSchedulerEndpoint: Attempted to request executors before the AM has registered!  
>WARN TableMapping: /etc/hadoop/conf/topology.table.file cannot be read.  
>java.io.FileNotFoundException: /etc/hadoop/conf/topology.table.file (No such file or directory)  
>...  
>WARN TableMapping: Failed to read topology table. /default-rack will be used for all nodes.  

and

>WARN Utils: Truncated the string representation of a plan since it was too large. This behavior can be adjusted by setting 'spark.debug.maxToStringFields' in SparkEnv.conf.  

which can be safely ignored.

### 4. Fit histograms

Histogram fitting uses local running or condor.

~~**Note:** the first time you run a fit it will compile the Root classes. Don't use `-j` option the first time you run fits.
It will try to compile the modules multiple times and throw errors. Instead, use single core to run one fit, then ctrl-c and use the `-j` option
(it won't compile again).~~

**(Update May 2021):** This extra compilation step shouldn't be needed anymore. Instead the fitting functions are compiled the first time the repository is sourced.

To run locally (with 16 threads):
```bash
./tnp_fitter.py fit muon generalTracks Z Run2017_UL configs/muon_pog_official_run2_Z_2017.json -j 16
```

To submit to condor:

First, do cd to your workspace in `condor_wrapper.sh` and give it execute permissions
```bash
$ sed -i "s~ReplaceMe_by_cdWorkdir~cd $PWD~" condor_wrapper.sh
$ chmod +x condor_wrapper.sh
```

Build the condor files and submit the jobs to the cluster:
```bash
$ ./tnp_fitter.py fit muon generalTracks Z Run2017_UL configs/muon_pog_official_run2_Z_2017.json --condor
$ condor_submit condor.sub
```

**Note:** CERN no longer allows condor jobs referencing eos mountpoints (see [here](https://batchdocs.web.cern.ch/troubleshooting/eos.html#no-eos-submission-allowed)). This means it's not possible to submit the jobs from the repository cloned on a SWAN instance, since SWAN projects are kept at `/eos/user/S/SOMEUSER/SWAN_projects`. So if you've cloned this repo into a SWAN area to look at the Jupyter notebooks, make sure to copy the repo elsewhere first before submitting condor jobs (e.g. afs user or work area).

**Note:** If you're running from the edge nodes, the condor submission step must still be done from lxplus (see above).

The histograms which are fit can be controlled with optional filters.
See documentation with:
```bash
./tnp_fitter.py fit -h
```
There is a simple automatic recovery processing that can be run
(in case of condor failures).
More advanced options (such as using statistical tests to evaluate the GOF)
are still being implemented.

```bash
./tnp_fitter.py fit muon generalTracks Z Run2017_UL configs/muon_pog_official_run2_Z_2017.json -j 16 --recover
```

**Note:** this recovery procedure simply looks for missing job output files in the fit directory, and runs the missing jobs again locally. If you are re-running your pipeline (e.g. after fixing some mistake in the config) using the same basedir, old fit jobs might still be located in that directory and so this simple recovery procedure won't detect the condor job failures. Additionally, you will end up with a mix of new and old fits in the folder. To avoid this, it's safer to erase the relevant folders before re-running (or choosing a new basedir).

### 5. Extract scale factors

Plots and scale factors can the be extracted with:

```bash
./tnp_fitter.py prepare muon generalTracks Z Run2017_UL configs/muon_pog_official_run2_Z_2017.json
```

Plots are saved into a plots/ directory and the json files into a efficiencies/ directory. Two json versions are now produced, one with the standard Muon POG format and another with the new xPOG recommended schemas. Feel free to use either.


**Note:** this is still a WIP.

## Utilities

### Pileup
The [make_pileup.py](make_pileup.py) script produced the pileup distribution in MC.
This part requires a CMSSW environment sourced.

To make the data pileup, copy the latest PileupHistogram from:
```bash
/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/{COLLISION_ERA}/{ENERGY}/PileUp/PileupHistogram-{...}.root
```
You should grab the `69200ub` version. If you wish to explore systematic uncertainties
in the choice of the minbias cross section, use the up (`66000ub`) and down (`72400ub`) histograms.

Alternatively, you can make it yourself with (e.g. Run2017):
```bash
lumimask=/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/ReReco/Cert_294927-306462_13TeV_EOY2017ReReco_Collisions17_JSON.txt
pileupjson=/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/PileUp/pileup_latest.txt
xsec=69200
maxBins=100
pileupCalc.py -i $lumimask --inputLumiJSON $pileupjson --calcMode true  --minBiasXsec $xsec --maxPileupBin $maxBins --numPileupBins $maxBins pileup/data/Run2017.root
```

### 6. Pre-fit distributions

Plot pre-fit distributions using local running or condor (add --condor_submit option). Some options are supported.

Plot Data/MC plots for a full era:

```
./tnp_fitter.py compare particle probe resonance era configs/muon_example.json --baseDir ./example
```

Produce plots for two sepecific suberas:

```
./tnp_fitter.py compare particle probe resonance era configs/muon_example.json --baseDir ./example --subera1 SubEra1 --subera2 SubEra2
```

If the two suberas came from different eras, just add the option "--era2":

```
./tnp_fitter.py compare particle probe resonance era1 configs/muon_example.json --baseDir ./example --subera1 SubEra1 --subera2 SubEra2 --era2 Era2
```
