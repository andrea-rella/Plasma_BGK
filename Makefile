# Compiler and flags
CXX      = g++-15

CXXFLAGS = -std=c++20 -Wall -Wextra

CPPFLAGS = -O3 \
    		-I/usr/local/Cellar/nlohmann-json/3.12.0/include \
            -I/usr/local/Cellar/eigen/3.4.0_1/include/eigen3 \
			-Iinclude


# Targets
EXEC = main.exe
SRC  = main.cpp
OBJ  = $(SRC:.cpp=.o)

# Default target
all: $(EXEC)

# Compile source files into object files
%.o: %.cpp
	$(CXX) -c $(CPPFLAGS) $(CXXFLAGS) $< -o $@

# Link object files into the final executable
$(EXEC): $(OBJ)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $^ $(LIBS) -o $@

# Run the program
run: $(EXEC)
	./$(EXEC)

# Clean up build artifacts
clean:
	$(RM) *.o

distclean: clean
	$(RM) *~ $(EXEC)


