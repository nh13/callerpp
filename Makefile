SRC_DIR=   src
BUILD_DIR= build
BIN_DIR=   bin
CPP=       g++
CPP_FLAGS= -std=c++11
OBJS=      $(BUILD_DIR)/caller.o
PROG=      callerpp
INCLUDES=
LIBS=


.PHONY: directories clean

.SUFFIXES:.cpp .o

$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp
	$(CPP) -c $(CPP_FLAGS) $(INCLUDES) $< -o $@

all:directories $(BIN_DIR)/$(PROG)

$(BIN_DIR)/$(PROG): $(OBJS)
	$(CPP) $(CPP_FLAGS) $(OBJS) -o $@ $(LIBS) -lspoa

$(BUILD_DIR) $(BIN_DIR):
	mkdir -p $@

directories: $(BUILD_DIR) $(BIN_DIR)

clean:
	rm -rv $(BUILD_DIR) $(BIN_DIR)
