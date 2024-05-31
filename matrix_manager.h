//############################################# About Author #########################################
// Created by: Alsamman M. Alsamman                                                                  #
// Emails: smahmoud [at] ageri.sci.eg or A.Alsamman [at] cgiar.org or SammanMohammed [at] gmail.com  #
// License: MIT License - https://opensource.org/licenses/MIT                                        #
// Disclaimer: The script comes with no warranty, use at your own risk                               #
// This script is not intended for commercial use                                                    #
//####################################################################################################


// create a 2D array of integers
int **createIntMatrix(int rows_size, int cols_size);
// create a 2D array of doubles
double **createDoubleMatrix(int rows_size, int cols_size);
// flip a matrix
int **flipMarkerMatrix(int **matrix, int numSamples, int numMarkers);
// print a 2D array of integers
void printIntMatrix(int **matrix, int rows_size, int cols_size);
////////////////////////////////////////// Memory Management //////////////////////////////////////////
// Free memory of a 2D array of integers
void freeIntMatrix(int **matrix, int rows_size);
// Free memory of a 2D array of doubles
void freeDoubleMatrix(double **matrix, int rows_size);
void binary2matrix(int **matrix, int marker_size, int sample_size, char *input_filename);
void matrix2binary(int **matrix, int marker_size, int sample_size, char *output_filename);
