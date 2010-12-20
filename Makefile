CC	= g++
CFLAGS	= -O2 -fopenmp
LIBS	= -lm -lpthread  

SRC=src/*.cpp
HDR=src/*.h

mstatx: $(SRC) $(HDR)
	@$(CC) $(CFLAGS) $(LIBS) -o mstatx $(SRC)

clean:
	@rm -f mstatx
