#INCDIR = -I$(ROOTSYS)/include -Iinc
INCDIR = -Iinc

LIBS = -lstdc++ -L/usr/local/lib -lm #$(shell root-config --evelibs) -lMLP -lMinuit -lTreePlayer -lTMVA -lTMVAGui -lXMLIO -lMLP -lm -lRooFit -lRooFitCore
CC = g++ -v

CFLAGS = -g -ansi -pedantic -W -Wall -Wextra -std=c++17 -m64
TARGET = CharmWidthDecay

SRC = $(wildcard src/*.cxx)

OBJECTS = $(subst src,obj,$(patsubst %.cxx,%.o,$(SRC)))

HEADERS = $(wildcard $(INCDIR)/*.h)

.PHONY: all clean

all : $(TARGET)
	@echo "DONE!"

obj/%.o: src/%.cxx $(HEADERS)
	$(CC) $(CFLAGS) $(INCDIR) -c $< -o $@

$(TARGET): $(OBJECTS) $(HEADERS)
	$(CC) $(CFLAGS) $(LIBS) $(OBJECTS) -o $@
obj:
	mkdir -p $@
clean:
	rm -rf ./obj
	rm -rf $(TARGET)
