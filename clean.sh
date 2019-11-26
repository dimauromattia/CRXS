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
rm -rf *.py *.pyc __pycache__  examples/example
popd


pushd $DIRNAME/py2/cpp
rm -rf __pycache__ *.pyc build  xs_tools.py xs_tools_wrap.cxx  _xs_tools*.so
rm -rf __pycache__ *.pyc build  xs_wrapper.py xs_wrapper_wrap.cxx  _xs_wrapper*.so
popd

pushd $DIRNAME/py2
rm -rf __pycache__ *.pyc XS_tools.py XS_wrapper.py
popd

pushd $DIRNAME/py3/cpp
rm -rf __pycache__ *.pyc build  xs_tools.py xs_tools_wrap.cxx  _xs_tools*.so
rm -rf __pycache__ *.pyc build  xs_wrapper.py xs_wrapper_wrap.cxx  _xs_wrapper*.so
popd

pushd $DIRNAME/py3
rm -rf __pycache__ *.pyc XS_tools.py XS_wrapper.py
popd
