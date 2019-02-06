#! /bin/bash

DIRNAME=`DIRNAME $0`
pushd $DIRNAME
FULLDIRNAME=`pwd`
popd

echo Linking package from $FULLDIRNAME!

#Linking internal:

ln -s $FULLDIRNAME/python/XS_wrapper.py  $FULLDIRNAME/XS_wrapper.py
ln -s $FULLDIRNAME/python/info.py        $FULLDIRNAME/info.py
ln -s $FULLDIRNAME/python/XS_tools.py    $FULLDIRNAME/XS_tools.py

echo " " >> $FULLDIRNAME/__init__.py


read -p "Specify your python3 sitepackage path:" spp
if [ -z $spp ]; then
  echo "Please specify a path."
  exit 1
fi

echo "Link to: $spp"
rm  $spp/CRXS
pushd $spp
ln -s $FULLDIRNAME CRXS
popd
