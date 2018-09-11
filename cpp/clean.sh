DIRNAME=`dirname $0`
echo Start cleaning from $DIRNAME!
rm $DIRNAME/include/*
rm $DIRNAME/lib/*
rm $DIRNAME/bin/*
rm -r $DIRNAME/build/*
