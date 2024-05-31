# compiler
CC = gcc

# compiler flags
CFLAGS = -Wall  -g
#-Werror
# linker flags
LDFLAGS = -lm

# # source files
SRCS = $(wildcard *.c)

# # object files
OBJS = $(SRCS:.c=.o)

# create the program

vcfCompressor: $(OBJS)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

main.o: main.c
	$(CC) $(CFLAGS) -c $< 

vcf_parser.o: vcf_parser.c
	$(CC) $(CFLAGS) -c $<

# run the program
run: vcfCompressor
	./vcfCompressor

# clean up
clean:
	rm -f vcfCompressor $(OBJS)

# valgrind
valgrind: vcfCompressor
	valgrind --leak-check=full ./vcfCompressor

test: vcfCompressor
	./vcfCompressor -c sampleData/sample.vcf
	./vcfCompressor -r -m sampleData/sample_markers.txt \
	-s sampleData/sample_samples.txt \
	-o sampleData/newvcf.vcf \
	-b sampleData/sample.bin