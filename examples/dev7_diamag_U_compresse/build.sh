rm -r build/
CC=mpicc CXX=mpicxx cmake . -B build
cmake --build build/
