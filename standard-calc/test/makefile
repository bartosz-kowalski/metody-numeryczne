CC = gcc
CFLAGS = -O2 

TARGET = test

all: $(TARGET)

test : test.o gradient.o 
	$(CC) -o $@ $^

test.o: test.c
	$(CC) $(CFLAGS) -c $< -o $@

gradient.o: ../standard-calc/gradient.c
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -f main test.o gradient.o

.PHONY: all clean
