
all:
	(mkdir -p ./cpp/build; cd ./cpp/build; cmake ..; make)
	(g++ -I./cpp/include -L./cpp/lib -lCRXS -o examples/example examples/example.cpp)
	(cd ./py2/cpp; swig -c++ -python xs_tools.i; python setup_xs_tools.py build_ext --inplace)
	(cd ./py2/cpp; swig -c++ -python xs_wrapper.i; python setup_xs_wrapper.py build_ext --inplace)
	(cd ./py3/cpp; swig -c++ -python -py3 xs_tools.i; python3 setup_xs_tools.py build_ext --inplace)
	(cd ./py3/cpp; swig -c++ -python -py3 xs_wrapper.i; python3 setup_xs_wrapper.py build_ext --inplace)

build:
	(mkdir -p ./cpp/build; cd ./cpp/build; cmake ..; make)

python:
	(cd ./py2/cpp; swig -c++ -python xs_tools.i; python setup_xs_tools.py build_ext --inplace)
	(cd ./py2/cpp; swig -c++ -python xs_wrapper.i; python setup_xs_wrapper.py build_ext --inplace)

python_three:
	(cd ./py3/cpp; swig -c++ -python -py3 xs_tools.i; python3 setup_xs_tools.py build_ext --inplace)
	(cd ./py3/cpp; swig -c++ -python -py3 xs_wrapper.i; python3 setup_xs_wrapper.py build_ext --inplace)

debug:
	(mkdir -p ./cpp/build; cd ./cpp/build; cmake -DCMAKE_BUILD_TYPE=Debug ..; make)

example:
	(g++ -I./cpp/include -L./cpp/lib -lCRXS -o examples/example examples/example.cpp)

clean:
	( bash clean.sh)
