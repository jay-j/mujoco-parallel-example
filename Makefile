COMMON=-O2 -Imujoco210/include -Lmujoco210/bin -g -std=c++11 -mavx -pthread -Wl,-rpath,'$$ORIGIN'
MUJOCO=mujoco210/bin/libmujoco210.so -lGL -lglew mujoco210/bin/libglfw.so.3 -llapacke
LIB=-Idependencies/
OBJLIB=mjgraphics.o

### Main Targets ###

all: main.bin

main.bin: $(OBJLIB) main.o
	g++ $(COMMON) $(OBJLIB) main.o $(MUJOCO) -o main.bin

### Libraries & Object Files ###
%.o: %.cpp 
	g++ $(COMMON) $(LIB) -c $< $(MUJOCO) -o $@

clean:
	rm *.o
	rm *.bin
