//############################################# About Author #########################################
// Created by: Alsamman M. Alsamman                                                                  #
// Emails: smahmoud [at] ageri.sci.eg or A.Alsamman [at] cgiar.org or SammanMohammed [at] gmail.com  #
// License: MIT License - https://opensource.org/licenses/MIT                                        #
// Disclaimer: The script comes with no warranty, use at your own risk                               #
// This script is not intended for commercial use                                                    #
//####################################################################################################
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
// our header files
#include "vcf_parser.h"


// reading vcf functions

void skipComments(FILE *vcfopen)
{
    while(1)
    {
        char c = fgetc(vcfopen);
        if(c == EOF)
        {
            printf("File is empty\n");
            exit(-1);
        }
        if(c == '#')
        {
            while(c != '\n')
            {
                c = fgetc(vcfopen);
            }
        }
        else
        {
            ungetc(c, vcfopen);
            break;
        }
    }
}
vcfinfo *readVCF(Parameters *params)
{
    ///////////////////////////
    vcfinfo *vcf = getBasicVCFInfo(params);
    // initialize the matrix in the vcfinfo
    vcf->matrix = createMarkerMatrix(vcf->marker_size, vcf->sample_size);
    vcf->markers_info = (marker_info *)malloc(vcf->marker_size * sizeof(marker_info));
    vcf->nchr = -1;
    vcf->sample_names = malloc(vcf->sample_size * sizeof(char *));
    for(int i = 0; i < vcf->sample_size; i++)
    {
        vcf->sample_names[i] = malloc(MAX_SAMPLE_NAME * sizeof(char));  
    }
    // reread the file and fill the matrix
    convertVCFToMatrix(vcf,params);
    // get the sample names
    getSampleNames(vcf, params);
    return vcf;
}
vcfinfo* initVCF(Parameters *params)
{
    vcfinfo *vcf = (vcfinfo *)malloc(sizeof(vcfinfo));
    vcf->nchr = -1;
    vcf->sample_size = params->sample_size;
    vcf->marker_size = params->marker_size;
    vcf->matrix = createMarkerMatrix(vcf->marker_size, vcf->sample_size);
    vcf->markers_info = (marker_info *)malloc(vcf->marker_size * sizeof(marker_info));
    vcf->sample_names = malloc(vcf->sample_size * sizeof(char *));
    for(int i = 0; i < vcf->sample_size; i++)
    {
        vcf->sample_names[i] = malloc(MAX_SAMPLE_NAME * sizeof(char));
    }
    return vcf;
}

// Memory management functions
void freeVCFInfo(vcfinfo *vcf)
{
    // free markers matrix
    freeIntMatrix(vcf->matrix, vcf->marker_size);
    // free sample names
    for(int i = 0; i < vcf->sample_size; i++)
    {
        free(vcf->sample_names[i]);
    }
    free(vcf->sample_names);
    // free the markers info
    free(vcf->markers_info);
    free(vcf);
}

vcfinfo *getBasicVCFInfo(Parameters *params)
{
    FILE *vcfopen = fopen(params->vcf_filename,"r");
    if(vcfopen == NULL)
    {
        fprintf(stderr, "Error: File not found\n");
        exit(1);
        return NULL;
    }
    skipComments(vcfopen);
    char c ;
    int col = 0;
    int sample_size = 0;
    int marker_size = 0;
    while((c = fgetc(vcfopen)) != EOF)
    {
        if(c == '\t' && sample_size == 0)
        {
            col++;
        }
        else if(c == '\n')
        {
            if(sample_size == 0)
                sample_size = col - 8;
            marker_size++;
        }
    }
    // add one Chromosome to test
    // malloc a vcf info
    vcfinfo *vcfptr = (vcfinfo *)malloc(sizeof(vcfinfo));
    vcfptr->nchr = -1;
    vcfptr->sample_size = sample_size;
    vcfptr->marker_size = marker_size;
    // printVCFInfo(vcf);
    fclose(vcfopen);
    return vcfptr;
}
int **createMarkerMatrix(int marker_size, int sample_size)
{
    int **matrix = (int **)malloc(marker_size * sizeof(int *));
    for(int i = 0; i < marker_size; i++)
    {
        matrix[i] = (int *)malloc(sample_size * sizeof(int));
    }
    // set the matrix to 9
    for(int i = 0; i < marker_size; i++)
    {
        for(int j = 0; j < sample_size; j++)
        {
            matrix[i][j] = 9;
        }
    }
    return matrix;
}
void addMarkerValue(int **matrix, int marker_size, int sample_size, int row, int col, char *value)
{
    if(row <= marker_size && col <= sample_size)
    {
        int value_int = -1;
        if(strstr(value, "0/0") != NULL)
            value_int = 0;
        else if(strstr(value, "0/1") != NULL)
            value_int = 1;
        else if(strstr(value, "1/0") != NULL)
            value_int = 1;
        else if(strstr(value, "1/1") != NULL)
            value_int = 2;
        else if(strstr(value, "./.") != NULL)
            value_int = MarkerMissing;
        else if(strstr(value, "0|0") != NULL)
            value_int = 0;
        else if(strstr(value, "0|1") != NULL)
            value_int = 1;
        else if(strstr(value, "1|0") != NULL)
            value_int = 1;
        else if(strstr(value, "1|1") != NULL)
            value_int = 2;
        else if(strstr(value, ".|.") != NULL)
            value_int = MarkerMissing;
        else
        {
            printf("Error: Unknown value %s\n", value);
            exit(-1);
            // ploy ploid markers will be treated as missing data
            value_int = MarkerMissing;
        }
        matrix[row][col] = value_int;
    }
    else
    {
        printf("Error: Out of range\n");
    }
}

void convertVCFToMatrix(vcfinfo *vcf, Parameters *params)
{
    FILE *vcfopen = fopen(params->vcf_filename, "r");
    if(vcfopen == NULL)
    {
        printf("File not found\n");
        exit(-1);
    }
    // read sample names
    skipComments(vcfopen);
    char c;
    int col = 0;
    int row = 0;
    char markerValue[100];
    int markerValueIndex = 0;

    // markers info
    int position = 0;
    char chr[MAX_CHR_NAME] = {0};
    char alleles[MAX_ALLELE+MAX_ALLELE+1] = {0};
    char rsid[MAX_MARKER_NAME] = {0};
    char allele1[MAX_ALLELE] = {0};
    char allele2[MAX_ALLELE] = {0};

    while((c = fgetc(vcfopen)) != EOF)
    {
        if(c == '\t' || c == '\n')
        {
            markerValue[markerValueIndex] = '\0';

            if (col == 0) // Chromosome column
            {
                strcpy(chr, markerValue);
            }
            else if (col == 1) // Position column
            {
                position = atoi(markerValue);
            }
            else if (col == 2) // rsID column
            {
                strcpy(rsid, markerValue);
            }
            else if (col == 3) // Allele1 column
            {
                strcpy(allele1, markerValue);
            }
            else if (col == 4) // Allele2 column
            {
                strcpy(allele2, markerValue);
            }
            else if (col >= 9) // Genotype information columns
            {
                addMarkerValue(vcf->matrix, vcf->marker_size, vcf->sample_size, row, col - 9, markerValue);
            }

            if (c == '\t')
            {
                col++;
            }
            else if (c == '\n')
            {
                // Store the marker information
                storeMarkerInfo(vcf, row, chr, position, rsid, allele1, allele2);
                // Reset marker information for the next row
                memset(chr, 0, MAX_CHR_NAME);
                memset(rsid, 0, MAX_MARKER_NAME);
                memset(alleles, 0, MAX_ALLELE);
                memset(allele1, 0, MAX_ALLELE);
                memset(allele2, 0, MAX_ALLELE);
                row++;
                col = 0;
            }
            memset(markerValue, 0, 100);
            markerValueIndex = 0;
        }
        else
        {
            markerValue[markerValueIndex] = c;
            markerValueIndex++;
        }
    }
    fclose(vcfopen);
}


void getSampleNames(vcfinfo *vcf, Parameters *params)
{
    FILE *vcfopen = fopen(params->vcf_filename, "r");
    if(vcfopen == NULL)
    {
        printf("File not found\n");
        exit(-1);
    }
    // loop until the #CHROM line, where its # then C then H then R then O then M
    while (1)
    {
        char c = fgetc(vcfopen);
        if (c == EOF)
        {
            printf("Error: No #CHROM line found\n");
            exit(-1);
        }
        if (c == '#')
        {
            char chrom[7] = {0};
            chrom[0] = c;
            int i;
            for (i = 1; i < 6; i++)
            {
                c = fgetc(vcfopen);
                // if c is not a letter, break
                if (c < 65 || c > 90|| c == EOF|| c == '\n'|| c == '\t')
                {
                    break;
                }
                chrom[i] = c;
            }
            chrom[i] = '\0';
            if (strcmp(chrom, "#CHROM") == 0)
            {
                break;
            }
        }
    }
    // read the sample names after column 9
    char c;
    int col = 0;
    int sample_index = 0;
    char sample_name[MAX_SAMPLE_NAME] = {0};
    while ((c = fgetc(vcfopen)) != EOF)
    {
        if (c == '\t' || c == '\n')
        {
            sample_name[sample_index] = '\0';
            if (col >= 9)
            {
                strcpy(vcf->sample_names[col - 9], sample_name);
            }
            if (c == '\t')
            {
                col++;
            }
            else if (c == '\n')
            {
                break;
            }
            memset(sample_name, 0, MAX_SAMPLE_NAME);
            sample_index = 0;
        }
        else
        {
            sample_name[sample_index] = c;
            sample_index++;
        }
    }
    fclose(vcfopen);


}

// This is a placeholder for the storeMarkerInfo function
void storeMarkerInfo(vcfinfo *vcf, int row, char *chr, int position, char *rsid, char *allele1, char *allele2)
{
    // Ensure the markers_info array is large enough
    if (row >= vcf->marker_size) {
        printf("Error: Row index %d exceeds marker_size %d\n", row, vcf->marker_size);
        exit(-1);
    }

    vcf->markers_info[row].marker_index = row;
    vcf->markers_info[row].position = position;
    strcpy(vcf->markers_info[row].chr, chr);
    strcpy(vcf->markers_info[row].rsid, rsid); // Assuming rsid is part of markers_info
    strcpy(vcf->markers_info[row].allele1, allele1);
    strcpy(vcf->markers_info[row].allele2, allele2);
}


void printMarkerMatrix(int **matrix, int marker_size, int sample_size)
{
    for(int i = 0; i < marker_size; i++)
    {
        for(int j = 0; j < sample_size; j++)
        {
            printf("%d", matrix[i][j]);
        }
        printf("\n");
    }
}
void printVCFInfo(vcfinfo *vcf)
{
    if(vcf->sample_size == 0)
        printf("Sample Size is not set\n");
    else
        printf("Sample Size: %d\n", vcf->sample_size);
    if(vcf->marker_size == 0)
        printf("Marker Size is not set\n");
    else
        printf("Marker Size: %d\n", vcf->marker_size);
}

int missingData(int *row, int size)
{
    int missing = 0;
    for (int i = 0; i < size; i++)
    {
        if (row[i] == MarkerMissing)
        {
            missing++;
        }
    }
    return missing;
}
void printMarkerInfo(vcfinfo *vcf, char *output_filename)
{
    FILE *output = fopen(output_filename, "w");
    if (output == NULL)
    {
        fprintf(stderr, "Error: File not found\n");
        exit(1);
    }
    fprintf(output, "Number of samples: %d\n", vcf->sample_size);
    fprintf(output, "Number of markers: %d\n", vcf->marker_size);
    marker_info *markers_info = vcf->markers_info;
    for (int i = 0; i < vcf->marker_size; i++)
    {
        fprintf(output, "%s %d %s %s %s\n", markers_info[i].chr, markers_info[i].position, markers_info[i].rsid, markers_info[i].allele1, markers_info[i].allele2);
    }
    fclose(output);
    printf("VCF info has been written to vcf_info.txt\n");
}

void readMarkerInfo(vcfinfo *vcf, char *input_filename)
{
    FILE *input = fopen(input_filename, "r");
    if (input == NULL)
    {
        fprintf(stderr, "Error: File not found\n");
        exit(1);
    }

    char line[1000];
    int i = 0;
    while (fgets(line, sizeof(line), input))
    {
        // skip the first 2 lines
        if (i == 0 || i==1)
        {
            i++;
            continue;
        }
        int position;
        char chr[100], rsid[100], allele1[100], allele2[100];
        sscanf(line, "%s %d %s %s %s\n", chr, &position, rsid, allele1, allele2);
        storeMarkerInfo(vcf, i-2, chr, position, rsid, allele1, allele2);
        i++;
    }
    fclose(input);

    printf("VCF info has been read from vcf_info.txt\n");
}
// using matrix and markers info to reconstruct the vcf file
void reconstructVCF(vcfinfo *vcf, char *output_filename)
{
    FILE *output = fopen(output_filename, "w");
    if (output == NULL)
    {
        fprintf(stderr, "Error: File not found\n");
        exit(1);
    }
    // print the header
    fprintf(output, "##fileformat=VCFv4.2\n");
    fprintf(output, "##source=VCFtools\n");
    fprintf(output, "##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of Samples With Data\">\n");
    fprintf(output, "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">\n");
    fprintf(output, "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency\">\n");
    fprintf(output, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n");
    fprintf(output, "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">\n");
    fprintf(output, "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">\n");
    fprintf(output, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT");
    for (int i = 0; i < vcf->sample_size; i++)
    {
        fprintf(output, "\t%s", vcf->sample_names[i]);
    }
    fprintf(output, "\n");
    // print the markers
    marker_info *markers_info = vcf->markers_info;
    for (int i = 0; i < vcf->marker_size; i++)
    {


        // fprintf(output, "%s\t%d\t%s\t%s\t%s\t.\t.\t.\tGT", markers_info[i].chr, markers_info[i].position, markers_info[i].rsid, allele1, allele2);
        fprintf(output, "%s\t%d\t%s\t%s\t%s\t.\t.\t.\tGT", markers_info[i].chr, markers_info[i].position, markers_info[i].rsid, markers_info[i].allele1, markers_info[i].allele2);

        for (int j = 0; j < vcf->sample_size; j++)
        {
            int value = vcf->matrix[i][j];
            switch (value)
            {
            case 0:
                fprintf(output, "\t0/0");
                break;
            case 1:
                fprintf(output, "\t0/1");
                break;
            case 2:
                fprintf(output, "\t1/1");
                break;
            case 9:
                fprintf(output, "\t./.");
                break;
            default:
                fprintf(stderr, "Invalid character encountered: %c\n", value);
                // write the invalid character to the screen
                fclose(output);
            }
        }
        fprintf(output, "\n");
    }
    fclose(output);
    printf("VCF file has been reconstructed to %s\n", output_filename);
}

void writeSampleNames(vcfinfo *vcf, char *filename)
{
    FILE *file = fopen(filename , "w");
    if (file == NULL)
    {
        printf("Error opening file!\n");
        exit(1);
    }
    for (int i = 0; i < vcf->sample_size; i++)
    {
        fprintf(file, "%s\n", vcf->sample_names[i]);
    }
    fclose(file);
}