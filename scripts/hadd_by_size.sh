if [ "$#" -ne 1 ]; then
	echo "Script requires exactly 1 argument: path containing the ROOT files to merge. Exiting..."
	exit 0
fi

sizeInMB=`du -ms $1 | awk '{print int($1)}'` # | sed 's|[^0-9\.]||g'` # extract total folder size in GB
numOutFiles=$(( $sizeInMB / 1000 )) # divide total size by 1GB to obtain total number of out files
numOutFiles=$(( $numOutFiles + 1 )) # add one to account for integer division truncation
numInFiles=`ls $1/*.root | wc -l` # total number of input ROOT files
if [ $numInFiles -eq 0 ]; then
	echo "No input ROOT files found. Exiting..."
	exit 0
fi
numInFilesPerOutFile=$((numInFiles / numOutFiles))

for i in `seq 0 $(($numOutFiles-1))`; do
	offset=$(($i * numInFilesPerOutFile))
	offset=$(($offset + 1))
	uuid=$(uuidgen | cut -f1 -d'-')
	numParallelProcs=$(( ($numInFiles - $i * $numInFilesPerOutFile) / 3 ))
	numParallelProcs=$(( $numParallelProcs < 17 ? $numParallelProcs : 16 ))
	if [ $numParallelProcs -eq 0 ]; then
		echo hadd haddOut_${i}_${uuid}.root `ls $1/*.root | tail -n +$offset | head -n $numInFilesPerOutFile`
		hadd haddOut_${i}_${uuid}.root `ls $1/*.root | tail -n +$offset | head -n $numInFilesPerOutFile`
	else
		hadd -j $numParallelProcs haddOut_${i}_${uuid}.root `ls $1/*.root | tail -n +$offset | head -n $numInFilesPerOutFile`
	fi
done
