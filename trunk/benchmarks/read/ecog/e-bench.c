#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <string.h>
#include "hdf5.h"
#include "e-bench.h"

int main(int argc, char* argv[]) {
    
    herr_t      status;
    hid_t       file_id, ECoGData_id, ECoGData_space, ECoGData_memspace, plist_id;
    hsize_t     i, j, ECoGData_dim[MAXDIM], my_count[MAXDIM], my_offset[MAXDIM], my_trials, total_trials;
    hsize_t*    read_idx;
    
    ECoGMeta*   metadata;
    
    char        query_labels[QUERY_SIZE][LABEL_LEN];
    char        filename[HNAME_MAX];

    int         query_num, cio=1;
    int         proc_num, my_rank;
    
    hsize_t     ECoGIndx_size;
    double*     ECoGIndx;
    double*     ECoGData_data;

    // Start MPI
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &proc_num);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    strcpy(filename, DATAFILE);

    // check arguments
    query_num = parse_argv(argc, argv, filename, query_labels, &cio);


    // Only root process open and read all index datasets
    if (my_rank == ROOTPROC) { 

        printf("%s\n", filename);
        // Allocate matadata storage 
        metadata     = (ECoGMeta*)malloc(sizeof(ECoGMeta));
        status       = root_get_metadata(filename, metadata);
 
        // Allocate memory for actual read index
        read_idx = (hsize_t*)malloc(metadata->EIndx_size* sizeof(hsize_t));
        total_trials = 0;
        // Find how much data we need to read
        for (i = 0; i < metadata->EIndx_dim[0]; i++) {
            for (j = 0; j < query_num; j++) {
                //printf("[%d] %s - %s\n", (int)(metadata->EIndx_data[i]) - 1, metadata->ELbls_data[ (int)(metadata->EIndx_data[i]) - 1 ], query_labels[j]);
                if ( strcmp( metadata->ELbls_data[ (int)(metadata->EIndx_data[i]) - 1 ], query_labels[j]) == 0 ) {
                    read_idx[total_trials++] = i;
                }
            }
        }
        ECoGIndx_size = metadata->ECoGIndx_size;

        //test(metadata);
        //get_indx_freq(metadata->EIndx_data, metadata->EIndx_size, metadata);

        printf("Total trials to be read %llu\n", total_trials);
    }
    
    // total_trials is the number of data blocks(each with 301*256 elements) that we need to read
    // read_idx has all the offsets

    MPI_Bcast(&total_trials, 1, MPI_LONG_LONG_INT, ROOTPROC, MPI_COMM_WORLD);

    // TODO: further decompose when proc_num is larger than needed trials.
    if(total_trials < proc_num) {
        if(my_rank == ROOTPROC)
            printf("Too few reads, aborting...\n");
        MPI_Finalize();
        exit(-1);
    }
    
    // Allocate memory for actual read index, all same size for now
    if(my_rank != ROOTPROC)
        read_idx = (hsize_t*)malloc(total_trials * sizeof(hsize_t));

    MPI_Bcast(read_idx, total_trials, MPI_LONG_LONG_INT, ROOTPROC, MPI_COMM_WORLD);
    
    // Broadcast index array
    MPI_Bcast(&ECoGIndx_size, 1, MPI_LONG_LONG_INT, ROOTPROC, MPI_COMM_WORLD);
    ECoGIndx = (double*)malloc(ECoGIndx_size * sizeof(double));
    //printf("ECoGIndx_size%d\n", ECoGIndx_size);
    if(my_rank == ROOTPROC)
        memcpy(ECoGIndx, metadata->ECoGIndx_data, ECoGIndx_size * sizeof(double));
    MPI_Bcast(ECoGIndx, ECoGIndx_size, MPI_DOUBLE, ROOTPROC, MPI_COMM_WORLD);
    
    // All processes open the data file.
    plist_id           = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL);

    // Open file and dataset
    file_id            = H5Fopen(filename, H5F_ACC_RDONLY, plist_id);
    ECoGData_id        = H5Dopen(file_id, DSET_NAME, H5P_DEFAULT);
    ECoGData_space     = H5Dget_space(ECoGData_id);

    // Get real data dimension
    H5Sget_simple_extent_dims(ECoGData_space, ECoGData_dim, NULL);

    if (my_rank != proc_num-1) {
        my_trials      = total_trials/proc_num;
    }
    else {
        my_trials      = total_trials/proc_num + total_trials%proc_num;
    }

    // Calculate total space needed 
    hsize_t tmp_idx, my_total_count[MAXDIM];
    my_total_count[0]  = 0; 
    my_total_count[1]  = ECoGData_dim[1]; 
    for (i = 0; i < my_trials; i++) {
        tmp_idx            = read_idx[(total_trials/proc_num)*my_rank+i] * 2;
        my_total_count[0] += ECoGIndx[tmp_idx+1] - ECoGIndx[tmp_idx] + 1; 
    }
    //printf("my_trials: %llu, my_total_count [0]:%llu, [1]:%llu\n", my_trials, my_total_count[0],my_total_count[1]);

    // Create memory space
    ECoGData_memspace  = H5Screate_simple(2, my_total_count, NULL);

    // File space 
    ECoGData_data      = (double*)malloc(my_total_count[0]*my_total_count[1]*sizeof(double));

    my_offset[1]   = 0;
    my_count[1]    = ECoGData_dim[1];

    for (i = 0; i < my_trials; i++) {
        tmp_idx        = read_idx[(total_trials/proc_num)*my_rank+i] * 2;
        my_offset[0]   = ECoGIndx[tmp_idx ] - 1; 
        my_count[0]    = ECoGIndx[tmp_idx+1] - ECoGIndx[tmp_idx] + 1; 

        if (i == 0)
            status     = H5Sselect_hyperslab(ECoGData_space, H5S_SELECT_SET, my_offset, NULL, my_count, NULL);
        else
            status     = H5Sselect_hyperslab(ECoGData_space, H5S_SELECT_OR, my_offset, NULL, my_count, NULL);
        //status = H5Sselect_hyperslab(ECoGData_space, H5S_SELECT_SET, my_offset, NULL, my_count, NULL);
        //status = H5Dread(ECoGData_id, H5T_IEEE_F64LE, ECoGData_memspace, ECoGData_space, H5P_DEFAULT, ECoGData_data+elem_trial*i);
        //printf("\n%d - Offsets: [0]:%llu [1]:%llu\n", my_rank, my_offset[0],my_offset[1]);
        //printf("%d - Counts:  [0]:%llu [1]:%llu\n", my_rank, my_count[0],my_count[1]);
    }

    double start_time, elapsed_time, all_time_max, all_time_min, all_time_avg;
    MPI_Barrier(MPI_COMM_WORLD);

    start_time = MPI_Wtime();

    // Are we doing collective read?
    if (cio) {
        hid_t plist2_id = H5Pcreate(H5P_DATASET_XFER);
        H5Pset_dxpl_mpio(plist2_id, H5FD_MPIO_COLLECTIVE);
        status = H5Dread(ECoGData_id, H5T_IEEE_F64LE, ECoGData_memspace, ECoGData_space, plist2_id, ECoGData_data);
        H5Pclose(plist2_id);
    }
    else {
        status = H5Dread(ECoGData_id, H5T_IEEE_F64LE, ECoGData_memspace, ECoGData_space, H5P_DEFAULT, ECoGData_data);
    }
    if (status < 0) {
        printf("Error reading! Exit...\n");
        exit(-1);
    }
    
    elapsed_time = MPI_Wtime() - start_time;

    MPI_Reduce(&elapsed_time, &all_time_max, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&elapsed_time, &all_time_min, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
    MPI_Reduce(&elapsed_time, &all_time_avg, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    all_time_avg /= proc_num;

    if(my_rank == ROOTPROC) {
         //printf("Total time: %f Min time: %f Avg time: %f Total data: %.1fM Agg Bandwidth: %f\n"
         //                , all_time_max, all_time_min, all_time_avg, total_trials*my_total_count[1]*sizeof(double)/1024.0/1024.0
         //                , my_total_count[0]*my_total_count[1]*sizeof(double)/1024.0/1024.0/all_time_max);
         printf("Total time: %f\t",all_time_max);
         printf("Total size: %lluM\n", total_trials*my_total_count[1]*sizeof(double)/1024.0/1024.0);
    }

    double my_sum = 0.0, all_sum;
    for (i = 0; i < my_total_count[0]*my_total_count[1]; i++) {
        // There are -nan exist in the datafile, so need to check
        if ( !isnan(ECoGData_data[i]) ) 
            my_sum += ECoGData_data[i]; 
    }
    
    MPI_Reduce(&my_sum, &all_sum, 1, MPI_DOUBLE, MPI_SUM, ROOTPROC, MPI_COMM_WORLD);

    if (my_rank == ROOTPROC) {
        
        // Sum used as checking correctness.
        printf("All sum = %f\n", all_sum);

        free(metadata->ECoGIndx_data);
        free(metadata->EIndx_data);
        for (i = 0; i < metadata->ELbls_dim[0]; i++) {
            free(metadata->ELbls_data[i]);
        }
        free(metadata->ELbls_data);
        free(metadata);
    }

    free(read_idx);
    free(ECoGData_data);
    free(ECoGIndx);

    // Close everything.
    status = H5Dclose(ECoGData_id);
    status = H5Sclose(ECoGData_space);
    status = H5Sclose(ECoGData_memspace);
    H5Pclose(plist_id);
    status = H5Fclose(file_id);

    MPI_Finalize();
    return 0;
}

herr_t root_get_metadata(char* filename, ECoGMeta* metadata) 
{
    int     i;
    hid_t   file_id;
    hid_t   ECoGIndx_id, EIndx_id, ELbls_id, ELbls_memtype;
    hid_t   ECoGIndx_space, EIndx_space, ELbls_space;
    herr_t  status;
     
    file_id = H5Fopen(filename, H5F_ACC_RDWR, H5P_DEFAULT);

    // Open datasets
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
    metadata->EIndx_size     = 1;
    // Calculate how much space we need
    for (i = 0; i < metadata->ECoGIndx_ndims; i++) 
        metadata->ECoGIndx_size *= metadata->ECoGIndx_dim[i];

    for (i = 0; i < metadata->EIndx_ndims; i++) 
        metadata->EIndx_size    *= metadata->EIndx_dim[i];

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

    // This is specially used for reading variable length strings.
    status = H5Dvlen_reclaim (ELbls_memtype, ELbls_space, H5P_DEFAULT, tmpELbls_data);
    free(tmpELbls_data);

    // Close
    status = H5Dclose(ECoGIndx_id);
    status = H5Dclose(EIndx_id);
    status = H5Dclose(ELbls_id);
    status = H5Sclose(ECoGIndx_space);
    status = H5Sclose(EIndx_space);
    status = H5Sclose(ELbls_space);
    status = H5Tclose(ELbls_memtype);
    status = H5Fclose(file_id);

    return status;
}

int parse_argv(int argc, char* argv[], char* filename, char query_labels[QUERY_SIZE][LABEL_LEN], int* cio)
{
    int i = 0;
    char c;
    while ((c = getopt (argc, argv, "f:q:c:")) != -1) {
        switch (c) 
        {   
        case 'f':
            strcpy(filename, optarg);
            break;
        case 'q':
            if (i >= QUERY_SIZE ) {
                printf("Increase QUERY_SIZE in .h file\nExit\n");
                exit(-1);
            }
            strcpy(query_labels[i++], optarg);
            break;
        case 'c':
            *cio = atoi(optarg);
            break;
        default:
            printf("Error option [%s]\n", optarg);
            printf("%s -f filename -c collective_io -q query_label\n", argv[0]);
            exit(-1);
        }   
    }  
    return i;
}

int test(ECoGMeta* metadata)
{
    printf("Event_ECoGIndx: [0]:%.0f [1]:%.0f\n", metadata->ECoGIndx_data[0], metadata->ECoGIndx_data[1]);
    printf("Event_ECoGIndx: [3510]:%.0f [3511]:%.0f\n", metadata->ECoGIndx_data[3510], metadata->ECoGIndx_data[3511]);
    printf("EIndx_data: [0]:%.0f [1755]:%.0f\n", metadata->EIndx_data[0], metadata->EIndx_data[1755]);
    printf("Labels: [0]:%s [50]:%s\n", metadata->ELbls_data[0], metadata->ELbls_data[50]);

    return 0;
}

int get_indx_freq(double* idx, int cnt, ECoGMeta* metadata)
{

    int i;
    int freq[LABEL_NUM];
    memset(freq, 0, sizeof(freq));
    for (i = 0; i < cnt; i++) {
        freq[ (int)(idx[i])-1 ]++;
    }
    for (i = 0; i < LABEL_NUM; i++) {
        printf("%s:%d \t", metadata->ELbls_data[i], freq[i]);
    }

    printf("\n");
}
