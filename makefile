OBJECTS = main.o
TARGET = main

# Compile
CXX = mpic++
CLIBS = $(which mpic++)
CFLAGS =
LDFLAGS =

# Run
CRUN = mpirun
NCPU = 4
CRFLAGS = -np $(NCPU)


#--------------------------------------------------
#			RULES
#--------------------------------------------------
all: $(TARGET)

$(TARGET): $(OBJECTS)
	$(CXX) $(OBJECTS) $(CFLAGS) -o $(TARGET) $(CLIBS)

clean:
	rm *.o

delete:
	rm *.o $(TARGET)

run:
	$(CRUN) $(CRFLAGS) $(TARGET)

# Rules to generate object files:
%.o: %.cpp
	$(CXX) $(CFLAGS) -c $<


