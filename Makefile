# Copyright (c) 2012 Guillaume Collet
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.

CC	= g++
CFLAGS	= -std=c++17 -O3 -Wall
LIBS	= -lm -lpthread

SRC=$(wildcard src/*.cpp)
SRC_NO_MAIN=$(filter-out src/main.cpp,$(SRC))
HDR=src/*.h

.PHONY: test clean

mstatx: $(SRC) $(HDR)
	$(CC) $(CFLAGS) $(LIBS) -o mstatx $(SRC)

# Every tests/test_XXX.cpp becomes its own standalone binary (its own main()),
# sharing tests/test_helpers.h. Add a new file here as new modules get covered.
TEST_BIN=tests/test_msa_scoring tests/test_jensen tests/test_kabat tests/test_wentropy tests/test_trident tests/test_gap tests/test_mvector tests/test_factory tests/test_options

test: $(TEST_BIN)

tests/test_msa_scoring: tests/test_msa_scoring.cpp tests/test_helpers.h $(SRC_NO_MAIN) $(HDR)
	$(CC) $(CFLAGS) $(LIBS) -I. -o tests/test_msa_scoring tests/test_msa_scoring.cpp $(SRC_NO_MAIN)
	./tests/test_msa_scoring

tests/test_jensen: tests/test_jensen.cpp tests/test_helpers.h $(SRC_NO_MAIN) $(HDR)
	$(CC) $(CFLAGS) $(LIBS) -I. -o tests/test_jensen tests/test_jensen.cpp $(SRC_NO_MAIN)
	./tests/test_jensen

tests/test_kabat: tests/test_kabat.cpp tests/test_helpers.h $(SRC_NO_MAIN) $(HDR)
	$(CC) $(CFLAGS) $(LIBS) -I. -o tests/test_kabat tests/test_kabat.cpp $(SRC_NO_MAIN)
	./tests/test_kabat

tests/test_wentropy: tests/test_wentropy.cpp tests/test_helpers.h $(SRC_NO_MAIN) $(HDR)
	$(CC) $(CFLAGS) $(LIBS) -I. -o tests/test_wentropy tests/test_wentropy.cpp $(SRC_NO_MAIN)
	./tests/test_wentropy

tests/test_trident: tests/test_trident.cpp tests/test_helpers.h $(SRC_NO_MAIN) $(HDR)
	$(CC) $(CFLAGS) $(LIBS) -I. -o tests/test_trident tests/test_trident.cpp $(SRC_NO_MAIN)
	./tests/test_trident

tests/test_gap: tests/test_gap.cpp tests/test_helpers.h $(SRC_NO_MAIN) $(HDR)
	$(CC) $(CFLAGS) $(LIBS) -I. -o tests/test_gap tests/test_gap.cpp $(SRC_NO_MAIN)
	./tests/test_gap

tests/test_mvector: tests/test_mvector.cpp tests/test_helpers.h $(SRC_NO_MAIN) $(HDR)
	$(CC) $(CFLAGS) $(LIBS) -I. -o tests/test_mvector tests/test_mvector.cpp $(SRC_NO_MAIN)
	./tests/test_mvector

tests/test_factory: tests/test_factory.cpp tests/test_helpers.h $(SRC_NO_MAIN) $(HDR)
	$(CC) $(CFLAGS) $(LIBS) -I. -o tests/test_factory tests/test_factory.cpp $(SRC_NO_MAIN)
	./tests/test_factory

tests/test_options: tests/test_options.cpp tests/test_helpers.h $(SRC_NO_MAIN) $(HDR)
	$(CC) $(CFLAGS) $(LIBS) -I. -o tests/test_options tests/test_options.cpp $(SRC_NO_MAIN)
	./tests/test_options

clean:
	rm -f mstatx tests/test_msa_scoring tests/test_jensen tests/test_kabat tests/test_wentropy tests/test_trident tests/test_gap tests/test_mvector tests/test_factory tests/test_options
