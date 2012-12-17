libpyplinkseqint: pyplinkseqint.cpp
	g++ -Wl,-no-undefined -Wl,--warn-once -fPIC -g -ggdb -rdynamic -I/usr/local/include/ -c pyplinkseqint.cpp -o pyplinkseqint.o
	g++ -Wl,-no-undefined -Wl,--warn-once -fPIC -g -ggdb -shared -rdynamic -L/usr/local/lib/ -o libpyplinkseqint.so pyplinkseqint.o -l plinkseq 