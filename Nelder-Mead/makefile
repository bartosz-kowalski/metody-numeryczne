CC = gcc
CFLAGS = -O2 

TARGET = simpleks

all: $(TARGET)

simpleks : main.o simpleks.o mean.o std.o
	$(CC) -o $@ $^

main.o: main.c
	$(CC) $(CFLAGS) -c $< -o $@

simpleks.o: Csimpleks.c
	$(CC) $(CFLAGS) -c $< -o $@

mean.o: ../standard-calc/mean.c
	$(CC) $(CFLAGS) -c $< -o $@

std.o: ../standard-calc/std.c
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -f main main.o simpleks.o mean.o std.o

.PHONY: all clean
