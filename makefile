# Compile
CXX       := mpic++
CXXFLAGS  := -pedantic-errors -Wall -Wextra -Werror
LDFLAGS   := -L/usr/lib -lstdc++ -lm
BUILD     := ./build
OBJ_DIR   := $(BUILD)/objects
EXE_DIR   := $(BUILD)/exec
TARGET    := program
INCLUDE   := -I include/
SRC       := $(wildcard src/*.cpp)
OBJECTS   := $(SRC:%.cpp=$(OBJ_DIR)/%.o)
INPUTFILE := ./IO/input
OUTPUTFILE:= ./IO/output
OBSFILE 	:= ./IO/obstacle

# Run
CRUN = mpirun
NCPU = 1
CRFLAGS = -np $(NCPU)


#--------------------------------------------------
#			RULES
#--------------------------------------------------
all: build $(EXE_DIR)/$(TARGET)

$(OBJ_DIR)/%.o: %.cpp
		@mkdir -p $(@D)
		$(CXX) $(CXXFLAGS) $(INCLUDE) -o $@ -c $<


$(EXE_DIR)/$(TARGET): $(OBJECTS)
		@mkdir -p $(@D)
		$(CXX) $(CXXFLAGS) $(INCLUDE) $(LDFLAGS) -o $(EXE_DIR)/$(TARGET) $(OBJECTS)

.PHONY: all build clean debug release


build:
		@mkdir -p $(EXE_DIR)
		@mkdir -p $(OBJ_DIR)

debug: CXXFLAGS += -DDEBUG -g
debug: all

release: CXXFLAGS += -O2
release: all

clean:
		-@rm -rvf $(OBJ_DIR)/*
		-@rm -rvf $(EXE_DIR)/*

run:
	$(CRUN) $(CRFLAGS) -v $(EXE_DIR)/$(TARGET) $(INPUTFILE) $(OUTPUTFILE) $(OBSFILE)
