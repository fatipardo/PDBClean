#!/bin/bash
# F.Poitevin - Stanford - June 2019
#
banner="----------\n OSF CLIENT\n ----------"
echo -e $banner
#
# ||||||||||||||||||||||||||
# ===== INITIALIZATION =====
# ||||||||||||||||||||||||||
#
# [source OSF]
if [ ! -f ~/.osfrc ]; then echo "Warning! ~/.osfrc was not found. This might prevent osf from working..."; else source ~/.osfrc; fi
# [inform the new user]
if [ $# -lt 3 ]; then  
  echo "Usage error: $0 <pull/push> <local directory> <remote directory> [optional: <update/static> <osf list filename (in static mode)> ]"
  echo ""
  echo "             !!! Mandatory arguments !!!"
  echo "                 - <pull/push>         : whether you want to upload or download files"
  echo "                 - <local directory>   : please give exact relative path"
  echo "                 - <remote directory>  : in pull mode, keyword-based search in list file. In push mode, it is ignored (write anything)."
  echo ""
  echo "             ~~~ Optional arguments ~~~"
  echo "                 - <update/static>     : whether to update the osf list file, or use an already existing one."
  echo "                 - <osf list filename> : name of local file where content of osf project will be listed"
  echo ""
  exit
fi
# [push or pull]
action=$1
# [check $localdir]
localdir=${2%/}
if [ ! -d $localdir ]; then echo "Usage error: $localdir does not exists."; exit; fi
# [check file list]
updatels=${4:-'no'}
if [ "$updatels" == 'yes' ]; then
  now=`date +"%Y%m%d-%H:%M"`
  osfls=${now}'_list.txt'
  echo "list of files on OSF will be stored in $osfls"
  osf ls > $osfls
else
  osfls=$5
  if [  ! -f $osfls ]; then
    echo "$osfls was not found. Updated list of files on OSF will be stored in $osfls"
    osf ls > $osfls
  else
    echo "list of files on OSF will be read from $osfls"
  fi
fi
echo "head $osfls : "
head $osfls
# [take care of remotedir now]
remotedir=${3%/}
if [$action == 'pull' ]; then
  if ! grep -q "$remotedir" $osfls; then
    echo "$remotedir was not found in OSF..."
    exit
  fi
else
  remotedir="/$(dirname "$localdir")"
fi
#
# |||||||||||||||||||||||
# ===== PUSH / PULL =====
# |||||||||||||||||||||||
#
if [ action == 'push' ]; then
  ### Scenario 1: directory does not exist in OSF
  if ! grep -q "$localdir" $osfls; then
    echo "$localdir is being created in OSF..."
    echo "> osf upload -r $localdir $remotedir"
    osf upload -r $localdir $remotedir
  ### Scenario 2: directory exists. Let's check it's complete
  else
    echo "uploading to $localdir in OSF..."
    for file in $localdir/*; do
      if ! grep -q "$file" $osfls; then
        echo ">>> $file ..."
        osf upload "$file" "/$file"; fi
    done
  fi
else
  echo "downloading content of $remotedir from OSF into current directory..."
  grep "$remotedir" $osfls > "tmp_${remotedir}.txt"
  while read -r line
  do
    if [ ! -f ${localdir}/${line##*/} ]; then
      echo ">>> $line ..."
      osf fetch $line ${localdir}/${line##*/}
    fi
  done < tmp_${remotedir}.txt
  rm -f tmp_${remotedir}.txt
fi
#
# |||||||||||||||||||
# ===== THE END =====
# |||||||||||||||||||
echo -e $banner

