CC=g++
CC+=-DDEBUG -g  
CFLAGS=-c -Wall -m64
LDFLAGS=-fPIC
DIR_SRC = ./src
SOURCES=Example.C $(wildcard $(DIR_SRC)/*.C)
OBJECTS=$(SOURCES:.C=.o)
EXECUTABLE=Example

CFLAGS += $(shell $(ROOTSYS)/bin/root-config --cflags)
LDFLAGS += $(shell $(ROOTSYS)/bin/root-config --libs) 

CFLAGS += -I./include/ -I$(ROOTSYS)/include/

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE):$(OBJECTS)
	$(CC) -o $@ $(OBJECTS) $(LDFLAGS)

.C.o:
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm -f *.o $(DIR_SRC)/*.o; rm $(EXECUTABLE)
