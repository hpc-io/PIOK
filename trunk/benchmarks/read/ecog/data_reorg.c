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

int data_reorg(char* filename, ECoGMeta* metadata)
{
    herr_t      status;
    hid_t       file_id, ECoGData_id, ECoGData_space, ECoGData_memspace, file_space;
    hid_t       new_ECoGData_id, new_ECoGData_space, new_ECoGData_memspace;
    hsize_t     i, j, k;
    hsize_t*    new_ECoGIndx = (hsize_t*)malloc(metadata->ECoGIndx_size*sizeof(hsize_t));

    hsize_t     ECoGData_dim[MAXDIM], file_sel[MAXDIM];
    hsize_t     my_count[MAXDIM], my_offset[MAXDIM];

    // Open file, dataset
    file_id              = H5Fopen(filename, H5F_ACC_RDWR, H5P_DEFAULT);
    ECoGData_id          = H5Dopen(file_id, DSET_NAME, H5P_DEFAULT);
    ECoGData_space       = H5Dget_space(ECoGData_id);


    H5Sget_simple_extent_dims(ECoGData_space, ECoGData_dim, NULL);

    my_offset[1]         = 0; 
    my_count[0]          = metadata->ECoGIndx_data[1] - metadata->ECoGIndx_data[0] + 1;
    my_count[1]          = ECoGData_dim[1]; 

    char new_dataname[128];
    char new_idxname[128];
    char new_chunkname[128];
    sprintf(new_dataname, "%s_hopt", DSET_NAME);
    sprintf(new_idxname, "%s_hidx", DSET_NAME);
    sprintf(new_chunkname, "%s_hchunk_size", DSET_NAME);
    printf("New dataset name: %s\n", new_dataname);

    // Create new dataset in file
    new_ECoGData_space   = H5Screate_simple(2, ECoGData_dim, NULL);
    new_ECoGData_id      = H5Dcreate(file_id, new_dataname, H5T_IEEE_F64LE, new_ECoGData_space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    

    // Allocate space for each trial
    hsize_t trial_size   = my_count[0] * my_count[1];
    double* ECoGData     = (double*)malloc(sizeof(double)*trial_size);

    printf("trial size:%llu EIndx_dim[0]=%llu\n", trial_size, metadata->EIndx_dim[0]);

    ECoGData_memspace    = H5Screate_simple(2, my_count, NULL);

    hsize_t write_offset = 0;
    // for all trial type
    hsize_t tmp;
    for (j = 0; j < metadata->ELbls_dim[0]; j++) {
   
        // for all trials
        for (i = 0; i < metadata->EIndx_dim[0]; i++) {
            tmp = (hsize_t)(metadata->EIndx_data[i]) - 1;
            if ( tmp == j ) {

                my_count[0]  = metadata->ECoGIndx_data[i*2+1] - metadata->ECoGIndx_data[i*2] + 1;
                my_offset[0] = (metadata->ECoGIndx_data[i*2] - 1);

                if (my_count[0] != 301) {
                    printf("%d: %f - %f = %llu\n", i, (metadata->ECoGIndx_data[i*2+1]), (metadata->ECoGIndx_data[i*2]), my_count[0]);
                }
                //printf("Read - Offsets: [0]:%llu [1]:%llu\n", my_offset[0],my_offset[1]);
                //printf("Read - Counts:  [0]:%llu [1]:%llu\n", my_count[0],my_count[1]);

                // read the data
                status       = H5Sselect_hyperslab(ECoGData_space, H5S_SELECT_SET, my_offset, NULL, my_count, NULL);
                status       = H5Dread(ECoGData_id, H5T_IEEE_F64LE, ECoGData_memspace, ECoGData_space, H5P_DEFAULT, ECoGData);
                //printf("[0]:%f,[1]:%f\n",ECoGData[0],ECoGData[1]);

                // update index
                // my_offset[0] is the offset of original layout
                // write_offset is the offset of new layout
                new_ECoGIndx[my_offset[0]/301 * 2    ]   = write_offset;
                new_ECoGIndx[my_offset[0]/301 * 2 + 1]   = 0;

                // write to new place
                my_offset[0] = write_offset;
                file_space   = H5Dget_space(new_ECoGData_id);
                status       = H5Sselect_hyperslab(file_space, H5S_SELECT_SET, my_offset, NULL, my_count, NULL);
                status       = H5Dwrite(new_ECoGData_id, H5T_IEEE_F64LE, ECoGData_memspace, file_space, H5P_DEFAULT, ECoGData);

                //printf("Write - Offsets: [0]:%llu [1]:%llu\n", my_offset[0],my_offset[1]);
                //printf("Write - Counts:  [0]:%llu [1]:%llu\n", my_count[0],my_count[1]);


                write_offset += my_count[0];
            }
            
        } // i

    } //j

    printf("k=%d\n",k);

    // Writing out reorganized layout index
    new_ECoGData_space   = H5Screate_simple(2, metadata->ECoGIndx_dim, NULL);
    new_ECoGData_id      = H5Dcreate(file_id, new_idxname, H5T_NATIVE_HSIZE, new_ECoGData_space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status               = H5Dwrite(new_ECoGData_id, H5T_NATIVE_HSIZE, H5S_ALL, H5S_ALL, H5P_DEFAULT, new_ECoGIndx);

    hid_t   chunk_space, chunk_id;
    hsize_t chunk_size[MAXDIM], chunk_size_dim[MAXDIM];
    chunk_size[0]        = 301;
    chunk_size[1]        = 256;
    chunk_size_dim[0]    = 2;
    chunk_size_dim[1]    = 1;
    
    chunk_space          = H5Screate_simple(2, chunk_size_dim, NULL);
    chunk_id             = H5Dcreate(file_id, new_chunkname, H5T_NATIVE_HSIZE, chunk_space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status               = H5Dwrite(chunk_id, H5T_NATIVE_HSIZE, H5S_ALL, H5S_ALL, H5P_DEFAULT, chunk_size);

    status = H5Dclose(chunk_id);
    status = H5Sclose(chunk_space);

    status = H5Dclose(new_ECoGData_id);
    status = H5Sclose(new_ECoGData_space);

    status = H5Dclose(ECoGData_id);
    status = H5Sclose(ECoGData_space);
    status = H5Sclose(ECoGData_memspace);
    status = H5Fclose(file_id);


}

int main(int argc, char* argv[])
{
    ECoGMeta* metadata = (ECoGMeta*)malloc(sizeof(ECoGMeta));

    char* fname;

    if(argc == 1)
        fname = "EC6_CV.h5";
    else
        fname = argv[1];
    read_metadata(fname, metadata);

    data_reorg(fname, metadata);

    return 0;
}

