#ifndef EBENCH_H
#define EBENCH_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <string.h>

#define DATAFILE    "./EC6_CV.h5"
#define DSET_NAME   "/Data/ECoG"
#define ROOTPROC    0
#define MAXDIM      3
#define LABEL_LEN   8
#define LABEL_NUM   53
#define QUERY_SIZE  16
#define HNAME_MAX   128

typedef struct ECoGMeta {
    // Dimension
    hsize_t ECoGIndx_ndims;
    hsize_t EIndx_ndims;
    hsize_t ELbls_ndims;   

    // Data size of each dimension
    hsize_t ECoGIndx_dim[MAXDIM];       // Offset range of each trial
    hsize_t EIndx_dim[MAXDIM];          // Type ID of each trial
    hsize_t ELbls_dim[MAXDIM];          // Label name of each type ID

    hsize_t ECoGIndx_size;
    hsize_t EIndx_size;

    // actual metadata
    double* ECoGIndx_data;
    double* EIndx_data;
    char**  ELbls_data;

} ECoGMeta;

herr_t root_get_metadata(char* filename, ECoGMeta* metadata);
int parse_argv(int argc, char* argv[], char* filename, char query_labels[QUERY_SIZE][LABEL_LEN], int* cio);
int test(ECoGMeta* metadata);


#endif
