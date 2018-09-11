#! /bin/bash
rm -rf html
rm -rf soft
mkdir soft
git clone https://github.com/korsmeier/CRXS.git soft/
cd soft/documentation/
doxygen Doxyfile
cp -r html ../..
cd ../..
rm -rf soft

