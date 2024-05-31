#include <stdio.h>
#include <stdlib.h>
#include "parameters.h"

void printParameters(Parameters *params)
{
    printf("Parameters:\n");
    printf("VCF file name: %s\n", params->vcf_filename);
    printf("Sample Size: %d\n", params->sample_size);
    printf("Marker Size: %d\n", params->marker_size);
}