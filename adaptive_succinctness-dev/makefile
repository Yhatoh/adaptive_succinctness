compile_normal_bv:
	g++ -I ~/include -L ~/lib src/main.cpp src/RunEncoderBitVector.cpp src/RunEncoderBitVector.hpp util/utils.hpp util/utils.cpp util/tunstallCoder.hpp util/tunstallCoder.cpp util/mysort.c.o util/nqsort.c.o util/huffman_coder.hpp util/huffman_coder.cpp -o fast_test -O9 -std=c++11 -lsdsl -ldivsufsort -ldivsufsort64

compile_normal_sd:
	g++ -I ~/include -L ~/lib src/main.cpp src/RunEncoderSDArray.cpp src/RunEncoderSDArray.hpp util/utils.hpp util/utils.cpp util/tunstallCoder.hpp util/tunstallCoder.cpp util/mysort.c.o util/nqsort.c.o util/huffman_coder.hpp util/huffman_coder.cpp -o fast_test -O9 -std=c++11 -lsdsl -ldivsufsort -ldivsufsort64

compile_suffix_bv:
	g++ -I ~/include -L ~/lib src/main.cpp src/RunEncoderBitVector.cpp src/RunEncoderBitVector.hpp util/utils.hpp util/utils.cpp util/tunstallCoder.hpp util/tunstallCoderSuffixArray.cpp util/mysort.c.o util/nqsort.c.o util/huffman_coder.hpp util/huffman_coder.cpp -o fast_test -O9 -std=c++11 -lsdsl -ldivsufsort -ldivsufsort64

compile_suffix_sd:
	g++ -I ~/include -L ~/lib src/main.cpp src/RunEncoderSDArray.cpp src/RunEncoderSDArray.hpp util/utils.hpp util/utils.cpp util/tunstallCoder.hpp util/tunstallCoderSuffixArray.cpp util/mysort.c.o util/nqsort.c.o util/huffman_coder.hpp util/huffman_coder.cpp -o fast_test -O9 -std=c++11 -lsdsl -ldivsufsort -ldivsufsort64

