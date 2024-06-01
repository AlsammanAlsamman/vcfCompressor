#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include "vcf_parser.h"
#include "parameters.h"
#include "matrix_manager.h"

void readSampleNames(vcfinfo *vcf, char *filename){
    FILE *file = fopen(filename, "r");
    if (file == NULL)
    {
        fprintf(stderr, "Error opening file %s\n", filename);
        exit(1);
    }
    for (int i = 0; i < vcf->sample_size; i++)
    {
        fscanf(file, "%s", vcf->sample_names[i]);
    }
    fclose(file);
}

// to get the number of samples and markers from the marker info file
void getNSNM(char *filename, int *num_samples, int *num_markers){
    FILE *file = fopen(filename, "r");
    if (file == NULL)
    {
        fprintf(stderr, "Error opening file %s\n", filename);
        exit(1);
    }
    fscanf(file, "Number of samples: %d\n", num_samples);
    fscanf(file, "Number of markers: %d\n", num_markers);
    fclose(file);
}

// Function to get the filename without extension
// Function to get the filename without the last .vcf extension
char* getFilenameWithoutExtension(const char* filepath) {
    // Duplicate the filepath to work with a modifiable copy
    char* result = strdup(filepath);
    if (result == NULL) {
        fprintf(stderr, "Error: Memory allocation failed for filename\n");
        exit(1);
    }

    // Find the last occurrence of ".vcf" in the result
    char* extension = strrchr(result, '.');
    if (extension != NULL && strcmp(extension, ".vcf") == 0) {
        *extension = '\0'; // Replace the '.' with null terminator to remove the extension
    }

    return result;
}


void printUsage(const char* programName) {
    printf("Usage:\n");
    printf("  %s -c <vcf_filename> -ns             Compress the VCF file\n", programName);
    printf("  %s -r -b <binary_file> -m <marker_file> -s <sample_file> -o <output_vcf>   Reconstruct the VCF file\n", programName);
    printf("  %s -h | --help                       Display this help message\n", programName);
}

int main(int argc, char *argv[]) {
    int opt;
    int compress_flag = 0, reconstruct_flag = 0;
    char *vcf_filename = NULL;
    char *binary_filename_in = NULL;
    char *marker_filename_in = NULL;
    char *sample_filename_in = NULL;
    char *output_filename_in = NULL;
    int num_samples = 0;
    int num_markers = 0;

    if (argc < 2) {
        printUsage(argv[0]);
        exit(1);
    }

    while ((opt = getopt(argc, argv, "hrc:b:m:s:o:")) != -1) {
        switch (opt) {
            case 'h':
                printUsage(argv[0]);
                exit(0);
            case 'c':
                compress_flag = 1;
                vcf_filename = argv[optind - 1]; // Use the next argument directly
                break;
            case 'r':
                reconstruct_flag = 1;
                vcf_filename = argv[optind - 1]; // Use the next argument directly
                break;
            case 'm':
                marker_filename_in = optarg;
                break;
            case 'b':
                binary_filename_in = optarg;
                break;
            case 's':
                sample_filename_in = optarg;
                break;
            case 'o':
                output_filename_in = optarg;
                break;
            default:
                printUsage(argv[0]);
                exit(1);
        }
    }

    if (compress_flag && reconstruct_flag) {
        fprintf(stderr, "Error: Cannot use both compress and reconstruct options together.\n");
        printUsage(argv[0]);
        exit(1);
    }

    if (compress_flag) {
        if (vcf_filename == NULL) {
            fprintf(stderr, "Error: VCF filename is required for compression.\n");
            printUsage(argv[0]);
            exit(1);
        }

        printf("Compressing VCF file: %s\n", vcf_filename);
        char *output_filename = getFilenameWithoutExtension(vcf_filename);
        char binary_filename[MAX_FILENAME_LENGTH];
        sprintf(binary_filename, "%s.bin", output_filename);
        char sample_filename[MAX_FILENAME_LENGTH];
        sprintf(sample_filename, "%s_samples.txt", output_filename);
        char marker_filename[MAX_FILENAME_LENGTH];
        sprintf(marker_filename, "%s_markers.txt", output_filename);

        Parameters params = {
            .vcf_filename = vcf_filename,
            .output_filename = output_filename
        };

        vcfinfo *vcf = readVCF(&params);
        matrix2binary(vcf->matrix, vcf->marker_size, vcf->sample_size, binary_filename);
        writeSampleNames(vcf, sample_filename);
        printMarkerInfo(vcf, marker_filename);
        freeVCFInfo(vcf);
        free(output_filename);

    } else if (reconstruct_flag) {
        if (vcf_filename == NULL || binary_filename_in == NULL || marker_filename_in == NULL || sample_filename_in == NULL || output_filename_in == NULL) {
            // which of the files is missing
            if (vcf_filename == NULL) {
                fprintf(stderr, "Error: VCF filename is required for reconstruction.\n");
            }
            if (binary_filename_in == NULL) {
                fprintf(stderr, "Error: Binary filename is required for reconstruction.\n");
            }
            if (marker_filename_in == NULL) {
                fprintf(stderr, "Error: Marker filename is required for reconstruction.\n");
            }
            if (sample_filename_in == NULL) {
                fprintf(stderr, "Error: Sample filename is required for reconstruction.\n");
            }
            if (output_filename_in == NULL) {
                fprintf(stderr, "Error: Output filename is required for reconstruction.\n");
            }

            fprintf(stderr, "Error: All file inputs and parameters are required for reconstruction.\n");
            printUsage(argv[0]);
            exit(1);
        }

        printf("Reconstructing VCF file from binary data: %s\n", vcf_filename);

        // Get the number of samples and markers from the marker info file
        getNSNM(marker_filename_in, &num_samples, &num_markers);

        Parameters params = {
            .vcf_filename = vcf_filename,
            .output_filename = output_filename_in,
            .marker_size = num_markers,
            .sample_size = num_samples
        };
        printf("Number of samples: %d\n", num_samples);
        printf("Number of markers: %d\n", num_markers);

        vcfinfo *new_vcf = initVCF(&params);
        readMarkerInfo(new_vcf, marker_filename_in);
        readSampleNames(new_vcf, sample_filename_in);
        binary2matrix(new_vcf->matrix, new_vcf->marker_size, new_vcf->sample_size, binary_filename_in);
        // // print the matrix
        reconstructVCF(new_vcf, output_filename_in);
        freeVCFInfo(new_vcf);
    } else {
        fprintf(stderr, "Error: You must specify either compress (-c) or reconstruct (-r) operation.\n");
        printUsage(argv[0]);
        exit(1);
    }
    return 0;
}