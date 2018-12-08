
all: build

build:
	(mkdir -p ./cpp/build; cd ./cpp/build; cmake ..; make)
	(cd ./cpp; swig -c++ -python -py3 xs_tools.i; python3 setup_xs_tools.py build_ext --inplace)
	(cd ./cpp; swig -c++ -python -py3 xs_wrapper.i; python3 setup_xs_wrapper.py build_ext --inplace)

cpp_only:
	(mkdir -p ./cpp/build; cd ./cpp/build; cmake ..; make)

debug:
	(mkdir -p ./cpp/build; cd ./cpp/build; cmake -DCMAKE_BUILD_TYPE=Debug ..; make)
	(cd ./cpp; swig -c++ -python xs_tools.i; python3 setup_xs_tools.py build_ext --inplace)
	(cd ./cpp; swig -c++ -python xs_wrapper.i; python3 setup_xs_wrapper.py build_ext --inplace)

example:
	(g++ -I./cpp/include -L./cpp/lib -lCRXS -o example example.cpp)

clean:
	(cd ./cpp; bash clean.sh)
	(cd ./cpp;    rm -rf __pycache__ build  xs_tools.py xs_tools_wrap.cxx  _xs_tools*.so)
	(cd ./cpp;    rm -rf __pycache__ build  xs_wrapper.py xs_wrapper_wrap.cxx  _xs_wrapper*.so)
	(cd ./python; rm -rf __pycache__)
	(cd ./kdd18;  rm -rf __pycache__)
	(rm example)
