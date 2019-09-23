CC = emcc
CXX = em++
CFLAGS= -Wall -O3
CXXFLAGS= -O3 --std=gnu++17 -s WASM=1 -s USE_BOOST_HEADERS=1 -I ../include/
VPATH=../src/wkt-nest/

OBJECTS=bbpack.o bbpack-geometry.o wktio.o

default: all
all:
	$(MAKE) clean
	$(MAKE) main
	$(MAKE) exports

test: check
check: all
	./test/all.sh

main: main.o $(OBJECTS)
	$(CXX) $(CXXFLAGS) -s EXIT_RUNTIME=1 -o wkt-nest.js $(OBJECTS)

exports: wasm-exports.o $(OBJECTS)
	$(CXX) $(CXXFLAGS) -s EXPORTED_FUNCTIONS="['_nest']" -s 'EXTRA_EXPORTED_RUNTIME_METHODS=["cwrap"]' -s MODULARIZE_INSTANCE=1 -s EXPORT_ES6=1 -o wkt-nest-exports.js wasm-exports.o $(OBJECTS)

.PHONY: clean
clean:
	rm -f *.o *.wasm *.js
