LCG_RELEASE=LCG_97apython3
LCG_ARCH=x86_64-centos7-gcc8
source /cvmfs/sft.cern.ch/lcg/views/${LCG_RELEASE}/${LCG_ARCH}-opt/setup.sh

if [[ "$HOSTNAME" == *"ithdp"* ]]; then
  # edge node
  echo "Sourcing hadoop edge node environment..."
  source hadoop-setconf.sh analytix
  echo "Done!"
elif [[ "$HOSTNAME" == *"lxplus"* ]]; then
  # lxplus
  echo "Sourcing lxplus environment..."
  source /cvmfs/sft.cern.ch/lcg/etc/hadoop-confext/hadoop-swan-setconf.sh analytix
  echo "Done!"
else
  echo "[Warning] Environment can only be lxplus or the CERN hadoop edge nodes. See README for more details"
  # still source lxplus env in case this is inside a condor node
  source /cvmfs/sft.cern.ch/lcg/etc/hadoop-confext/hadoop-swan-setconf.sh analytix
fi

# Compile the Roofit fitting function if it doesn't exist yet
if [ ! -f RooCMSShape_cc.so ]; then
    echo ""
    echo "Did not detect a RooCMSShape shared object file. Compiling with ACLiC... (should be needed only once)"
    root -l -b -q -e '.L RooCMSShape.cc+'
    echo "Done!"
fi

if [ ! -f RooErfExp_cc.so ]; then
    echo ""
    echo "Did not detect a RooErfExp shared object file. Compiling with ACLiC... (should be needed only once)"
    root -l -b -q -e '.L RooErfExp.cc+'
    echo "Done!"
fi
