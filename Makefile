SRC_DIR=   src
BUILD_DIR= build
BIN_DIR=   bin
TEST_DIR=  tests
CXX=       g++
CXX_FLAGS= -std=c++14
OBJS=      $(BUILD_DIR)/caller.o
PROG=      callerpp
INCLUDES=
LIBS=

# Test configuration
TEST_OBJS=     $(BUILD_DIR)/test_main.o $(BUILD_DIR)/test_caller.o
TEST_PROG=     test_callerpp
TEST_LIBS=     -lgtest -lpthread


.PHONY: all directories clean test

.SUFFIXES:.cpp .o

$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp
	$(CXX) -c $(CXX_FLAGS) $(INCLUDES) $< -o $@

$(BUILD_DIR)/%.o: $(TEST_DIR)/%.cpp
	$(CXX) -c $(CXX_FLAGS) $(INCLUDES) $< -o $@

all:directories $(BIN_DIR)/$(PROG)

$(BIN_DIR)/$(PROG): $(OBJS)
	$(CXX) $(CXX_FLAGS) $(OBJS) -o $@ $(LIBS) -lspoa

$(BIN_DIR)/$(TEST_PROG): $(TEST_OBJS)
	$(CXX) $(CXX_FLAGS) $(TEST_OBJS) -o $@ $(TEST_LIBS)

$(BUILD_DIR) $(BIN_DIR):
	mkdir -p $@

directories: $(BUILD_DIR) $(BIN_DIR)

test: directories $(BIN_DIR)/$(TEST_PROG)
	./$(BIN_DIR)/$(TEST_PROG)

clean:
	rm -rv $(BUILD_DIR) $(BIN_DIR)
