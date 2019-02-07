#! /bin/bash

DIRNAME=`DIRNAME $0`
pushd $DIRNAME
FULLDIRNAME=`pwd`
popd

echo Linking package from $FULLDIRNAME!

read -p "Please specify which python version you want to link (2 or 3):" version
read -p "Specify your python sitepackage path:" spp
if [ -z $spp ]; then
  echo "Please specify a path."
  exit 1
fi

linkdir=$FULLDIRNAME
if [ "$version" == "2" ]; then
   linkdir=$linkdir/py2
elif [ "$version" == "3" ]; then
   linkdir=$linkdir/py3
else
  echo "Please enter either 2 or 3 as version."
  exit 1
fi


echo "Link to: $spp"
rm  $spp/CRXS
pushd $spp
ln -s $linkdir CRXS
popd
