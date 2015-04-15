all: cms

CXX := g++

INCLUDE = .
CFLAGS = -O3 -Wall

OBJ_CPU = data.o main.o load.o util.o rmsd.o main.o

%.o : %.C
	$(CXX) -I$(INCLUDE) $(CFLAGS) -c $<

cms : $(OBJ_CPU)
	$(CXX) $^ -o $@
	@mkdir -p ../bin
	@echo "mv cms ../bin"
	@mv cms ../bin

clean :
	rm ./*.o
	rm ../bin/cms
