compile_normal_bv:
	g++ -I ~/include -L ~/lib src/main.cpp src/RunEncoderAccess.cpp src/RunEncoderAccess.hpp util/utils.hpp util/utils.cpp util/tunstallCoder.hpp util/tunstallCoder.cpp util/mysort.c.o util/nqsort.c.o util/huffman_coder.hpp util/huffman_coder.cpp -o fast_test -O3 -std=c++17 -lsdsl -ldivsufsort -ldivsufsort64

compile_normal_sd:
	g++ -I ~/include -L ~/lib src/main.cpp src/RunEncoderAccess.cpp src/RunEncoderAccess.hpp src/RunEncoderSelect.cpp src/RunEncoderSelect.hpp util/utils.hpp util/utils.cpp util/tunstallCoder.hpp util/tunstallCoder.cpp util/mysort.c.o util/nqsort.c.o util/huffman_coder.hpp util/huffman_coder.cpp -o fast_test -O3 -std=c++17 -lsdsl -ldivsufsort -ldivsufsort64

compile_suffix_bv:
	g++ -I ~/include -L ~/lib src/main.cpp src/RunEncoderAccess.cpp src/RunEncoderAccess.hpp util/utils.hpp util/utils.cpp util/tunstallCoder.hpp util/tunstallCoderSuffixArray.cpp util/mysort.c.o util/nqsort.c.o util/huffman_coder.hpp util/huffman_coder.cpp -o fast_test -O3 -std=c++17 -lsdsl -ldivsufsort -ldivsufsort64

compile_suffix_sd:
	g++ -I ~/include -L ~/lib src/main.cpp src/RunEncoderSelect.cpp src/RunEncoderSelect.hpp util/utils.hpp util/utils.cpp util/tunstallCoder.hpp util/tunstallCoderSuffixArray.cpp util/mysort.c.o util/nqsort.c.o util/huffman_coder.hpp util/huffman_coder.cpp -o fast_test -O3 -std=c++17 -lsdsl -ldivsufsort -ldivsufsort64

compile_test_normal_bv:
	g++ -I ~/include -L ~/lib src/test_access.cpp src/RunEncoderAccess.cpp src/RunEncoderAccess.hpp util/utils.hpp util/utils.cpp util/tunstallCoder.hpp util/tunstallCoder.cpp util/mysort.c.o util/nqsort.c.o util/huffman_coder.hpp util/huffman_coder.cpp -o test_bitvector -O3 -std=c++17 -lsdsl -ldivsufsort -ldivsufsort64

compile_test_normal_sd:
	g++ -pg -I ~/include -L ~/lib src/test_select.cpp src/RunEncoderAccess.cpp src/RunEncoderAccess.hpp src/RunEncoderSelect.cpp src/RunEncoderSelect.hpp util/utils.hpp util/utils.cpp util/tunstallCoder.hpp util/tunstallCoder.cpp util/mysort.c.o util/nqsort.c.o util/huffman_coder.hpp util/huffman_coder.cpp -o fast_testlast -O3 -std=c++17 -lsdsl -ldivsufsort -ldivsufsort64

compile_test_normal_pasta:
	g++-8 -I ~/include -I ~/adaptive_succinctness/bit_vector/include -I ~/adaptive_succinctness/bit_vector/extlib/utils/include -I ~/mnt/d/Universidad/adaptive_succinctness/bit_vector/extlib/tlx/ -L ~/lib src/main.cpp src/RunEncoderPasta.cpp src/RunEncoderPasta.hpp util/utils.hpp util/utils.cpp util/tunstallCoder.hpp util/tunstallCoder.cpp util/mysort.c.o util/nqsort.c.o util/huffman_coder.hpp util/huffman_coder.cpp -o test_pasta -O3 -lsdsl -ldivsufsort -ldivsufsort64 -std=c++2a -Wdeprecated-declarations

compile_test_suffix_sd:
	g++ -I ~/include -L ~/lib src/test_select.cpp src/RunEncoderSelect.cpp src/RunEncoderSelect.hpp util/utils.hpp util/utils.cpp util/tunstallCoder.hpp util/tunstallCoderSuffixArray.cpp util/mysort.c.o util/nqsort.c.o util/huffman_coder.hpp util/huffman_coder.cpp -o fast_test -O3 -std=c++17 -lsdsl -ldivsufsort -ldivsufsort64
