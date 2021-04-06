CC=gcc

CFLAGS=-Wall
OFLAGS=-O3

all: bin/genysis bin/split bin/split2

bin/genysis: src/genysis.c src/load.h src/detection.h src/popcount.h src/rdtsc.h
	$(CC) $(CFLAGS) $(OFLAGS) $< -o $@

bin/split: src/split.c
	$(CC) $(CFLAGS) $(OFLAGS) $< -o $@

bin/split2: src/split2.c
	$(CC) $(CFLAGS) $(OFLAGS) $< -o $@

clean:
	rm -Rf *~ src/*~ bin/genysis bin/split bin/split2

clean_file:
	rm 6* 7* L* MN* MT* MW* N* description.txt
