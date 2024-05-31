#ifndef PARAMETERS_H
#define PARAMETERS_H

typedef struct  {
    char *vcf_filename; // vcf file name
    char *output_filename; // output file name
    int sample_size; // number of samples
    int marker_size; // number of markers
} Parameters;

void printParameters(Parameters *params);
#endif