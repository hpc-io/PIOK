#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <string.h>
#include "hdf5.h"
#include "e-bench.h"

int read_metadata(char* filename, ECoGMeta* metadata)
{
    int     i;
    hid_t   file_id;
    hid_t   ECoGIndx_id, EIndx_id, ELbls_id, ELbls_memtype;
    hid_t   ECoGIndx_space, EIndx_space, ELbls_space;
    herr_t      status;

    file_id = H5Fopen(filename, H5F_ACC_RDWR, H5P_DEFAULT);

    // Open datasets, assuming the file is already open and file_id is correct
    ECoGIndx_id = H5Dopen(file_id, "/Descriptors/Event_ECoGIndx", H5P_DEFAULT);
    EIndx_id    = H5Dopen(file_id, "/Descriptors/Event_EIndx", H5P_DEFAULT);
    ELbls_id    = H5Dopen(file_id, "/Descriptors/Event_ELbls", H5P_DEFAULT);

    // Get space
    ECoGIndx_space = H5Dget_space(ECoGIndx_id);
    EIndx_space    = H5Dget_space(EIndx_id   );
    ELbls_space    = H5Dget_space(ELbls_id   );

    // Get data size and dim info
    metadata->ECoGIndx_ndims = H5Sget_simple_extent_dims(ECoGIndx_space, metadata->ECoGIndx_dim, NULL);
    metadata->EIndx_ndims    = H5Sget_simple_extent_dims(EIndx_space,    metadata->EIndx_dim, NULL);
    metadata->ELbls_ndims    = H5Sget_simple_extent_dims(ELbls_space,    metadata->ELbls_dim, NULL);

    metadata->ECoGIndx_size  = 1;
    // Calculate how much space we need
    for (i = 0; i < metadata->ECoGIndx_ndims; i++)
        metadata->ECoGIndx_size *= metadata->ECoGIndx_dim[i];

    metadata->EIndx_size    = metadata->EIndx_dim[0];

    // Allocate memory
    metadata->ECoGIndx_data = (double*)malloc(metadata->ECoGIndx_size*sizeof(double));
    metadata->EIndx_data    = (double*)malloc(metadata->EIndx_size*sizeof(double));

    char** tmpELbls_data    = (char**)malloc(metadata->ELbls_dim[0]*sizeof(char*));

    metadata->ELbls_data    = (char**)malloc(metadata->ELbls_dim[0]*sizeof(char*));
    for (i = 0; i < metadata->ELbls_dim[0]; i++) {
        metadata->ELbls_data[i] = (char*)malloc(LABEL_LEN*sizeof(char));
    }

    // Read all metadata
    status = H5Dread(ECoGIndx_id, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, metadata->ECoGIndx_data);
    status = H5Dread(EIndx_id   , H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, metadata->EIndx_data);

    // Special case for reading label strings.
    ELbls_memtype = H5Tcopy(H5T_C_S1);
    status = H5Tset_size(ELbls_memtype, H5T_VARIABLE);
    status = H5Dread(ELbls_id, ELbls_memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, tmpELbls_data);
    for (i = 0; i < metadata->ELbls_dim[0]; i++) {
        strcpy(metadata->ELbls_data[i], tmpELbls_data[i]);
    }


    // Close
    status = H5Dvlen_reclaim (ELbls_memtype, ELbls_space, H5P_DEFAULT, tmpELbls_data);
    free(tmpELbls_data);
    status = H5Dclose(ECoGIndx_id);
    status = H5Dclose(EIndx_id);
    status = H5Dclose(ELbls_id);
    status = H5Sclose(ECoGIndx_space);
    status = H5Sclose(EIndx_space);
    status = H5Sclose(ELbls_space);
    status = H5Tclose(ELbls_memtype);
    status = H5Fclose(file_id);
}

int data_expend(char* filename, ECoGMeta* metadata, hsize_t increase_mul)
{
    herr_t      status;
    hid_t       file_id, ECoGData_id, ECoGData_space, ECoGData_memspace, file_space;
    char        new_filename[HNAME_MAX];
    hid_t       new_file_id, new_ECoGData_id, new_ECoGData_g, new_ECoGData_space, new_ECoGData_memspace;
    hsize_t     i, j, k;


    hsize_t     ECoGData_dim[MAXDIM], file_sel[MAXDIM];
    hsize_t     new_ECoGData_dim[MAXDIM], new_ECoGIndx_dim[MAXDIM], new_EIndx_dim[MAXDIM];
    hsize_t     my_count[MAXDIM], my_offset[MAXDIM];

    // New filename
    sprintf(new_filename,"%s_%dG", filename, increase_mul);
    printf("%s\n", new_filename);

    // Open old file, dataset
    file_id              = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
    ECoGData_id          = H5Dopen(file_id, DSET_NAME, H5P_DEFAULT);
    ECoGData_space       = H5Dget_space(ECoGData_id);

    // Get dim
    H5Sget_simple_extent_dims(ECoGData_space, ECoGData_dim, NULL);

    // Space for all ECoG data
    double* ECoGData     = (double*)malloc(ECoGData_dim[0]*ECoGData_dim[1]*sizeof(double));

    // Read all ECoG_data of old file
    status               = H5Dread(ECoGData_id, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, ECoGData);
        
    //printf("ECoGData: [0]:%f [1]:%f [2]:%f\n",ECoGData[0],ECoGData[1],ECoGData[2]);

    // Init new dim
    new_ECoGData_dim[0]  = ECoGData_dim[0]*increase_mul;
    new_ECoGData_dim[1]  = ECoGData_dim[1];

    new_ECoGIndx_dim[0]  = metadata->ECoGIndx_dim[0]*increase_mul;
    new_ECoGIndx_dim[1]  = metadata->ECoGIndx_dim[1];

    new_EIndx_dim[0]     = metadata->EIndx_dim[0]*increase_mul;
    new_EIndx_dim[1]     = metadata->EIndx_dim[1];

    // Open new file, dataset
    new_file_id          = H5Fcreate(new_filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    new_ECoGData_g       = H5Gcreate(new_file_id, "/Data", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    new_ECoGData_space   = H5Screate_simple(2, new_ECoGData_dim, NULL);
    new_ECoGData_id      = H5Dcreate(new_ECoGData_g, "/Data/ECoG", H5T_IEEE_F64LE, new_ECoGData_space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    file_space           = H5Dget_space(new_ECoGData_id);


    // Space for modified ECoGIndx
    double*     new_ECoGIndx = (double*)malloc(metadata->ECoGIndx_size*increase_mul*sizeof(double));
    double*     new_EIndx    = (double*)malloc(metadata->EIndx_size*increase_mul*sizeof(double));
    if(new_ECoGIndx == NULL || new_EIndx == NULL) {
        printf("Malloc failed\n");
        exit(-1);
    }

    my_offset[1]         = 0; 

    hsize_t new_ECoGIndx_off, new_EIndx_off;

    for (j = 0; j < increase_mul; j++) {

        my_offset[0]     = j*ECoGData_dim[0]; 
        new_ECoGIndx_off = j*metadata->ECoGIndx_size;
        new_EIndx_off    = j*metadata->EIndx_size;
        printf("Start to append at location %llu, ECoGIndx_off=%llu\n", my_offset[0], new_ECoGIndx_off);

        status           = H5Sselect_hyperslab(file_space, H5S_SELECT_SET, my_offset, NULL, ECoGData_dim, NULL);
        status           = H5Dwrite(new_ECoGData_id, H5T_IEEE_F64LE, ECoGData_space, file_space, H5P_DEFAULT, ECoGData);

        // ECoGIndx
        for (i = 0; i < metadata->ECoGIndx_size; i++) {
            new_ECoGIndx[i+new_ECoGIndx_off] = metadata->ECoGIndx_data[i] + my_offset[0];
        }

        for (i = 0; i < metadata->EIndx_size; i++) {
            new_EIndx[i+new_EIndx_off]    = metadata->EIndx_data[i];
        }

    }


    new_ECoGData_g       = H5Gcreate(new_file_id, "/Descriptors", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    new_ECoGData_space   = H5Screate_simple(2, new_ECoGIndx_dim, NULL);
    new_ECoGData_id      = H5Dcreate(new_ECoGData_g, "/Descriptors/Event_ECoGIndx", H5T_IEEE_F64LE, new_ECoGData_space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status               = H5Dwrite(new_ECoGData_id, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, new_ECoGIndx);

    new_ECoGData_space   = H5Screate_simple(2, new_EIndx_dim, NULL);
    new_ECoGData_id      = H5Dcreate(new_ECoGData_g, "/Descriptors/Event_EIndx", H5T_IEEE_F64LE, new_ECoGData_space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status               = H5Dwrite(new_ECoGData_id, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, new_EIndx);

    hid_t       filetype = H5Tcopy (H5T_C_S1);
    status               = H5Tset_size (filetype, H5T_VARIABLE);
    hid_t memtype        = H5Tcopy (H5T_C_S1);
    status               = H5Tset_size (memtype, H5T_VARIABLE);

    hsize_t dims[1];
    dims[0] = metadata->ELbls_dim[0];
    hid_t space          = H5Screate_simple (1, dims, NULL);
    hid_t dset           = H5Dcreate (new_file_id, "/Descriptors/Event_ELbls", filetype, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status               = H5Dwrite (dset, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, metadata->ELbls_data);


    status = H5Dclose(new_ECoGData_id);
    status = H5Sclose(new_ECoGData_space);
    status = H5Gclose(new_ECoGData_g);

    status = H5Dclose(ECoGData_id);
    status = H5Sclose(ECoGData_space);
    
    status = H5Fclose(file_id);
    status = H5Fclose(new_file_id);


}

int main(int argc, char* argv[])
{
    ECoGMeta* metadata = (ECoGMeta*)malloc(sizeof(ECoGMeta));

    char* fname;
    int   mul;

    if(argc < 3) {
        fname =  "EC6_CV.h5";
        mul   = 5;
    }
    else if (argc == 3) {
        fname = argv[1];
        mul   = atoi(argv[2]); 
    }
    read_metadata(fname, metadata);

    data_expend(fname, metadata, mul);

    return 0;
}

