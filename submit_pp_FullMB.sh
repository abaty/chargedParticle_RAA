if [ $# -ne 0 ]
then
  echo "Usage: ./psort.sh <trackqual> <file-list> <tag> <nmin> <nmax> <pttrigmin> <pttrigmax> <ptassmin> <ptassmax>"
  exit 1
fi

now="Rpp_$(date +"%Y_%m_%d__%H_%M_%S")"
njobs=1800

mkdir $now
cp fileLists/testFileList.txt $now 
cp fileLists/pp_MB_VanderbiltAllPDs.txt $now/fileList.txt 
cp -r TrkCorr_* $now
cp -r chi2Corrector $now
tar -cvzf $now/Corrections.tar.gz TrkCorr_* chi2Corrector
cp run.sh $now

cat run.condor | sed "s@log_flag@$now@g" | sed "s@dir_flag@$PWD/$now@g" | sed "s@user_flag@$USER@g" |  sed "s@arglist@$njobs 1@g" | sed "s@transfer_filelist@run.exe@g" | sed "s@njobs@$njobs@g" > $now/run.condor

cat run.sh | sed "s@/temporaryStorage/@/temporaryStorage/ppFromVanderbilt/@g" > $now/run.sh

NAME="countTracks.C"
g++ $NAME $(root-config --cflags --libs) -Wall -O2 -o "run.exe"
cp run.exe $now
rm run.exe
echo finished compilation
echo
cat $now/run.condor
echo 
echo condor_submit $now/run.condor
echo
# condor_submit $now/run.condor

