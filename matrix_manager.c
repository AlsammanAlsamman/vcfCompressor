//############################################# About Author #########################################
// Created by: Alsamman M. Alsamman                                                                  #
// Emails: smahmoud [at] ageri.sci.eg or A.Alsamman [at] cgiar.org or SammanMohammed [at] gmail.com  #
// License: MIT License - https://opensource.org/licenses/MIT                                        #
// Disclaimer: The script comes with no warranty, use at your own risk                               #
// This script is not intended for commercial use                                                    #
//####################################################################################################
#include <stdio.h>
#include <stdlib.h>
#include "matrix_manager.h"

int **createIntMatrix(int rows_size, int cols_size)
{
    int **matrix = (int **)malloc(sizeof(int *) * rows_size);
    size_t i;
    for(i = 0; i < rows_size; i++)
    {
        matrix[i] = (int *)malloc(sizeof(int) * cols_size);
    }
    return matrix;
}

double **createDoubleMatrix(int rows_size, int cols_size)
{
    double **matrix = (double **)malloc(sizeof(double *) * rows_size);
    size_t i;
    for(i = 0; i < rows_size; i++)
    {
        matrix[i] = (double *)malloc(sizeof(double) * cols_size);
    }
    return matrix;
}


int **flipMarkerMatrix(int **matrix, int numSamples, int numMarkers)
{
    // create a temporary matrix
    int **tempMatrix = (int **)malloc(sizeof(int *) * numSamples);
    for (int i = 0; i < numSamples; i++)
    {
        tempMatrix[i] = (int *)malloc(sizeof(int) * numMarkers);
    }

    // flip the matrix
    for (int i = 0; i < numSamples; i++)
    {
        for (int j = 0; j < numMarkers; j++)
        {
            tempMatrix[i][j] = matrix[j][i];
        }
    }

    return tempMatrix;
}

////////////////////////////////////////// Memory Management //////////////////////////////////////////
void freeIntMatrix(int **matrix, int rows_size)
{
    size_t i;
    for(i = 0; i < rows_size; i++)
    {
        free(matrix[i]);
    }
    free(matrix);
}

void freeDoubleMatrix(double **matrix, int rows_size)
{
    size_t i;
    for(i = 0; i < rows_size; i++)
    {
        free(matrix[i]);
    }
    free(matrix);
}

void printIntMatrix(int **matrix, int rows_size, int cols_size)
{
    for(int i = 0; i < rows_size; i++)
    {
        for(int j = 0; j < cols_size; j++)
        {
            printf("%d", matrix[i][j]);
        }
        printf("\n");
    }
}

void matrix2binary(int **matrix, int marker_size, int sample_size, char *output_filename)
{
    FILE *output = fopen(output_filename, "wb");
    if (output == NULL)
    {
        fprintf(stderr, "Error: File not found\n");
        exit(1);
    }
    // the values will be written as follows:
    // 0: 00
    // 1: 01
    // 2: 10
    // 9: 11
    // bit by bit
    unsigned char byte = 0;
    int bit_count = 0;
    for (int i = 0; i < marker_size; i++)
    {
        for (int j = 0; j < sample_size; j++)
        {
            int value = matrix[i][j];
            switch (value) {
            case 0:
                byte = (byte << 1) | 0;
                bit_count++;
                if (bit_count == 8) {
                    fwrite(&byte, sizeof(byte), 1, output);
                    byte = 0;
                    bit_count = 0;
                }
                byte = (byte << 1) | 0;
                bit_count++;
                if (bit_count == 8) {
                    fwrite(&byte, sizeof(byte), 1, output);
                    byte = 0;
                    bit_count = 0;
                }
                break;
            case 1:
                byte = (byte << 1) | 0;
                bit_count++;
                if (bit_count == 8) {
                    fwrite(&byte, sizeof(byte), 1, output);
                    byte = 0;
                    bit_count = 0;
                }
                byte = (byte << 1) | 1;
                bit_count++;
                if (bit_count == 8) {
                    fwrite(&byte, sizeof(byte), 1, output);
                    byte = 0;
                    bit_count = 0;
                }
                break;
            case 2:
                byte = (byte << 1) | 1;
                bit_count++;
                if (bit_count == 8) {
                    fwrite(&byte, sizeof(byte), 1, output);
                    byte = 0;
                    bit_count = 0;
                }
                byte = (byte << 1) | 0;
                bit_count++;
                if (bit_count == 8) {
                    fwrite(&byte, sizeof(byte), 1, output);
                    byte = 0;
                    bit_count = 0;
                }
                break;
            case 9:
                byte = (byte << 1) | 1;
                bit_count++;
                if (bit_count == 8) {
                    fwrite(&byte, sizeof(byte), 1, output);
                    byte = 0;
                    bit_count = 0;
                }
                byte = (byte << 1) | 1;
                bit_count++;
                if (bit_count == 8) {
                    fwrite(&byte, sizeof(byte), 1, output);
                    byte = 0;
                    bit_count = 0;
                }
                break;
            default:
                fprintf(stderr, "Invalid character encountered: %c\n", value);
                // write the invalid character to the screen
                fclose(output);
        }
        }
    }
    // if there are any remaining bits
    if (bit_count > 0)
    {
        byte = byte << (8 - bit_count);
        fwrite(&byte, sizeof(byte), 1, output);
    }


    fclose(output);
    printf("Binary matrix has been written to %s\n", output_filename);
}


void binary2matrix(int **matrix, int marker_size, int sample_size, char *input_filename)
{
    FILE *input = fopen(input_filename, "rb");
    if (input == NULL)
    {
        fprintf(stderr, "Error: File not found\n");
        exit(1);
    }

    // Calculate total number of bits needed
    int total_bits = marker_size * sample_size * 2;
    // Calculate total number of bytes needed
    int total_bytes = (total_bits + 7) / 8; // Add 7 to round up to the nearest byte

    // Allocate memory to read the file contents
    unsigned char *buffer = (unsigned char *)malloc(total_bytes);
    if (buffer == NULL)
    {
        fprintf(stderr, "Error: Memory allocation failed\n");
        fclose(input);
        exit(1);
    }

    // Read the exact number of bytes
    size_t bytes_read = fread(buffer, 1, total_bytes, input);
    if (bytes_read != total_bytes)
    {
        fprintf(stderr, "Error: Could not read the expected number of bytes\n");
        free(buffer);
        fclose(input);
        exit(1);
    }
    fclose(input);

    int bit_count = 0;
    int value = 0;
    int i = 0, j = 0;
    int bits_processed = 0;

    for (int byte_index = 0; byte_index < total_bytes; byte_index++)
    {
        unsigned char byte = buffer[byte_index];
        for (int k = 7; k >= 0; k--)
        {
            int bit = (byte >> k) & 1;

            value = (value << 1) | bit;
            bit_count++;
            bits_processed++;

            if (bit_count == 2)
            {
                if (i >= marker_size || j >= sample_size)
                {
                    fprintf(stderr, "Error: Exceeded matrix dimensions\n");
                    free(buffer);
                    exit(1);
                }

                switch (value)
                {
                case 0b00:
                    matrix[i][j] = 0;
                    break;
                case 0b01:
                    matrix[i][j] = 1;
                    break;
                case 0b10:
                    matrix[i][j] = 2;
                    break;
                case 0b11:
                    matrix[i][j] = 9;
                    break;
                default:
                    fprintf(stderr, "Error: Invalid binary value\n");
                    free(buffer);
                    exit(1);
                }
                value = 0;
                bit_count = 0;
                j++;
                if (j == sample_size)
                {
                    j = 0;
                    i++;
                }
            }

            // Stop if the total number of bits processed equals the total bits needed
            if (bits_processed == total_bits)
            {
                free(buffer);
                printf("Matrix has been reconstructed from %s\n", input_filename);
                return;
            }
        }
    }

    free(buffer);
    printf("Matrix has been reconstructed from %s\n", input_filename);
}