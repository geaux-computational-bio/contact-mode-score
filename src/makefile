CXX := g++


OBJ_CPU := coords.o data.o main.o


%.o: %.C
	$(CXX) -c $<

main : $(OBJ_CPU)
	$(CXX) $^ -o $@

clean :
	rm ./*.o
	rm ./main
