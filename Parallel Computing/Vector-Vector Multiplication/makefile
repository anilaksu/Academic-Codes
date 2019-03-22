# Declaration of variables
CC = mpicxx
CC_FLAGS = -w
LDFLAGS = -lcblas -llapacke 

# File names
EXEC = MatrixVector
SOURCES = MatrixVector.cpp
OBJECTS = $(SOURCES:.cpp=.o)
 
# Main target
$(EXEC): $(OBJECTS)
	$(CC) $(OBJECTS) -o $(EXEC) $(LDFLAGS) 
 
# To obtain object files
%.o: %.cpp
	$(CC) -c $(CC_FLAGS) $< -o $@
 
# To remove generated files
clean:
	rm -f $(EXEC) $(OBJECTS)