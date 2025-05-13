CC = gcc
CFLAGS = -O2 

TARGET = newton

all: $(TARGET)

newton : main.o gradient.o newton.o
	$(CC) -o $@ $^

main.o: main.c
	$(CC) $(CFLAGS) -c $< -o $@

gradient.o: ../standard-calc/gradient.c
	$(CC) $(CFLAGS) -c $< -o $@

newton.o: newton.c
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -f main main.o newton.o gradient.o

.PHONY: all clean
