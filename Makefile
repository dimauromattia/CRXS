
all: build

build:
	(mkdir -p ./cpp/build; cd ./cpp/build; cmake ..; make)
	(cd ./cpp; swig -c++ -python -py3 xs_tools.i; python3 setup_xs_tools.py build_ext --inplace)
	(cd ./cpp; swig -c++ -python -py3 xs_wrapper.i; python3 setup_xs_wrapper.py build_ext --inplace)
	(g++ -I./cpp/include -L./cpp/lib -lCRXS -o examples/example examples/example.cpp)

cpp_only:
	(mkdir -p ./cpp/build; cd ./cpp/build; cmake ..; make)

debug:
	(mkdir -p ./cpp/build; cd ./cpp/build; cmake -DCMAKE_BUILD_TYPE=Debug ..; make)
	(cd ./cpp; swig -c++ -python xs_tools.i; python3 setup_xs_tools.py build_ext --inplace)
	(cd ./cpp; swig -c++ -python xs_wrapper.i; python3 setup_xs_wrapper.py build_ext --inplace)
	(g++ -I./cpp/include -L./cpp/lib -lCRXS -o examples/example examples/example.cpp)

example:
	(g++ -I./cpp/include -L./cpp/lib -lCRXS -o examples/example examples/example.cpp)

clean:
	( bash clean.sh)
