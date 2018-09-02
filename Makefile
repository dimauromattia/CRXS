
all: build

build:
	(cd ./cpp; swig -c++ -python clike.i; python3 setup.py build_ext --inplace)

clean:
	(cd ./cpp;    rm -rf __pycache__ build  clike.py clike_wrap.cxx  _clike*.so)
	(cd ./python; rm -rf __pycache__)
	(cd ./kdd18;  rm -rf __pycache__)
