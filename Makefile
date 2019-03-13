
all:
	(mkdir -p ./cpp/build; cd ./cpp/build; cmake ..; make)
	(g++ -I./cpp/include -L./cpp/lib -lCRXS -o examples/example examples/example.cpp)
	(cd ./py2/cpp; swig -c++ -python xs_tools.i; python setup_xs_tools.py build_ext --inplace)
	(cd ./py2/cpp; swig -c++ -python xs_wrapper.i; python setup_xs_wrapper.py build_ext --inplace)
	(cp python/XS_tools.py.in py2/XS_tools.py; sed -i "s/@INPORT_crxs@/cpp.xs_tools/" py2/XS_tools.py)
	(cp python/XS_wrapper.py.in py2/XS_wrapper.py; sed -i "s/@INPORT_crxs@/cpp.xs_wrapper/" py2/XS_wrapper.py)
	(cd ./py3/cpp; swig -c++ -python -py3 xs_tools.i; python3 setup_xs_tools.py build_ext --inplace)
	(cd ./py3/cpp; swig -c++ -python -py3 xs_wrapper.i; python3 setup_xs_wrapper.py build_ext --inplace)
	(cp python/XS_tools.py.in py3/XS_tools.py; sed -i "s/@INPORT_crxs@/CRXS.cpp.xs_tools/" py3/XS_tools.py)
	(cp python/XS_wrapper.py.in py3/XS_wrapper.py; sed -i "s/@INPORT_crxs@/CRXS.cpp.xs_wrapper/" py3/XS_wrapper.py)

build:
	(mkdir -p ./cpp/build; cd ./cpp/build; cmake ..; make)

python:
	(cd ./py2/cpp; swig -c++ -python xs_tools.i; python setup_xs_tools.py build_ext --inplace)
	(cd ./py2/cpp; swig -c++ -python xs_wrapper.i; python setup_xs_wrapper.py build_ext --inplace)
	(cp python/XS_tools.py.in py2/XS_tools.py; sed -i "s/@INPORT_crxs@/cpp.xs_tools/" py2/XS_tools.py)
	(cp python/XS_wrapper.py.in py2/XS_wrapper.py; sed -i "s/@INPORT_crxs@/cpp.xs_wrapper/" py2/XS_wrapper.py)

python_three:
	(cd ./py3/cpp; swig -c++ -python -py3 xs_tools.i; python3 setup_xs_tools.py build_ext --inplace)
	(cd ./py3/cpp; swig -c++ -python -py3 xs_wrapper.i; python3 setup_xs_wrapper.py build_ext --inplace)
	(cp python/XS_tools.py.in py3/XS_tools.py; sed -i "s/@INPORT_crxs@/CRXS.cpp.xs_tools/" py3/XS_tools.py)
	(cp python/XS_wrapper.py.in py3/XS_wrapper.py; sed -i "s/@INPORT_crxs@/CRXS.cpp.xs_wrapper/" py3/XS_wrapper.py)

debug:
	(mkdir -p ./cpp/build; cd ./cpp/build; cmake -DCMAKE_BUILD_TYPE=Debug ..; make)

example:
	(g++ -I./cpp/include -L./cpp/lib -lCRXS -o examples/example examples/example.cpp)

clean:
	( bash clean.sh)
