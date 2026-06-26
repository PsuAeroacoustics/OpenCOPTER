# Simple Makefile for OpenCOPTER C example

CC = gcc
#CFLAGS = -Wall -Wextra -std=c11 -g3
CFLAGS = -Wall -Wextra -std=c11 -O3
INCLUDES = -I../../include
LIBS = -L../../ -lopencopter -Xlinker -rpath $(PWD)/../../ -lm

TARGET = example
SOURCE = example.c

$(TARGET): $(SOURCE)
	$(CC) $(CFLAGS) $(INCLUDES) -o $(TARGET) $(SOURCE) $(LIBS)

clean:
	rm -f $(TARGET)

.PHONY: clean