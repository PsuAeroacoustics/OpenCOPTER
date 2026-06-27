# Simple Makefile for OpenCOPTER C example

CC = g++
#CFLAGS = -Wall -Wextra -std=c11 -g3
CFLAGS = -Wall -Wextra -std=c++23 -O3
INCLUDES = -I../../include
LIBS = -L../../ -lopencopter -Xlinker -rpath $(PWD)/../../ -lm

TARGET = example
SOURCE = example.cpp

$(TARGET): $(SOURCE)
	$(CC) $(CFLAGS) $(INCLUDES) -o $(TARGET) $(SOURCE) $(LIBS)

clean:
	rm -f $(TARGET)

.PHONY: clean