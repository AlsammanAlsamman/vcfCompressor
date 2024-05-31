//############################################# About Author #########################################
// Created by: Alsamman M. Alsamman                                                                  #
// Emails: smahmoud [at] ageri.sci.eg or A.Alsamman [at] cgiar.org or SammanMohammed [at] gmail.com  #
// License: MIT License - https://opensource.org/licenses/MIT                                        #
// Disclaimer: The script comes with no warranty, use at your own risk                               #
// This script is not intended for commercial use                                                    #
//####################################################################################################

#include <stdio.h>
#include "globals.h"
#include "matrix_manager.h"
#include "parameters.h"

#ifndef VCF_PARSER_H
#define VCF_PARSER_H
// Global variables
#define MAX_CHR_NAME 100
#define MAX_ALLELE 100
#define MAX_SAMPLE_NAME 100
#define MAX_MARKER_NAME 100
#define MarkerMissing 9
#define MAX_FILENAME_LENGTH 200

typedef struct{
    int marker_index;
    int position;
    char rsid[MAX_MARKER_NAME];
    char chr[MAX_CHR_NAME];
    char allele1[MAX_ALLELE];
    char allele2[MAX_ALLELE];
} marker_info;

// Data structures
typedef struct {
    int nchr;
    int sample_size;
    int marker_size;
    char **sample_names;
    marker_info *markers_info;
    int **matrix; // row: sample, col: marker
} vcfinfo;

#endif

// Function prototypes
// initialization functions
vcfinfo* initVCF(Parameters *params);
// printig functions
void printMarkerMatrix(int **matrix, int marker_size, int sample_size);
// reading vcf functions
void skipComments(FILE *vcfopen);
vcfinfo *getBasicVCFInfo(Parameters *params);
// void addMarkerInfo(vcfinfo *vcf, int marker_index, int position, char *chr, char *alleles);
int **createMarkerMatrix(int marker_size, int sample_size);
void addMarkerValue(int **matrix, int marker_size, int sample_size, int row, int col, char *value);
void convertVCFToMatrix(vcfinfo *vcf, Parameters *params);
vcfinfo *readVCF(Parameters *params);
void printVCFInfo(vcfinfo *vcf);
void storeMarkerInfo(vcfinfo *vcf, int row, char *chr, int position, char *rsid, char *allele1, char *allele2);
void freeVCFInfo(vcfinfo *vcf);
void getSampleNames(vcfinfo *vcf, Parameters *params);
void reconstructVCF(vcfinfo *vcf, char *output_filename);
void readMarkerInfo(vcfinfo *vcf, char *input_filename);
void printMarkerInfo(vcfinfo *vcf, char *output_filename);
void writeSampleNames(vcfinfo *vcf, char *filename);