DIRNAME=`dirname $0`
echo Start cleaning from $DIRNAME!

pushd $DIRNAME/cpp
bash clean.sh
rm -rf __pycache__ *.pyc build  xs_tools.py xs_tools_wrap.cxx  _xs_tools*.so
rm -rf __pycache__ *.pyc build  xs_wrapper.py xs_wrapper_wrap.cxx  _xs_wrapper*.so
popd
pushd $DIRNAME/python
rm -rf __pycache__ *.pyc
popd
pushd $DIRNAME
rm -rf *.py *.pyc __pycache__  example
popd
