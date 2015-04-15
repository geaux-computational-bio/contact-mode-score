CXX := g++


OBJ_CPU := data.o main.o load.o util.o rmsd.o
CFLAGS := -O3 -Wall

%.o : %.C
	$(CXX) $(CFLAGS) -c $<

cms : $(OBJ_CPU)
	$(CXX) $(CFLAGS) $^ -o $@

clean :
	rm ./*.o
	rm ./cms
