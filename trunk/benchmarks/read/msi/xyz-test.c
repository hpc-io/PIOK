//#include "sds-vol.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include "hdf5.h"

#define NAME_MAX 255
#define DIM 3

char filename[NAME_MAX];
char dataset_name[NAME_MAX];
char cb_buffer_size[NAME_MAX];
char cb_nodes[NAME_MAX];
char global_query[NAME_MAX];

int ht_query_parser(const char *query, hsize_t *x_min, hsize_t *x_max, hsize_t *y_min, hsize_t *y_max, hsize_t *z_min, hsize_t *z_max) {

    char tmp[512];
    char* token;
    // "[100:200,200:400,200:400]"

    sscanf(query, "[%[^]]]", tmp);

    token = strtok(tmp, ",");
    if(*token == ':') {
        *x_min  = 0;
        *x_max  = -1; 
    }   
    else {
        sscanf(token,"%llu:%llu",x_min, x_max);
    }   

    token = strtok(NULL, ",");
    if(*token == ':') {
        *y_min  = 0;
        *y_max  = -1; 
    }   
    else {
        sscanf(token,"%llu:%llu",y_min, y_max);
    }   

    token = strtok(NULL, ",");
    if(*token == ':') {
        *z_min  = 0;
        *z_max  = -1; 
    }   
    else {
        sscanf(token,"%llu:%llu",z_min, z_max);
    }   


    return 0;
}


int main(int argc, char *argv[]) {

    int        i, j, k, c, mpi_size, mpi_rank, cio = 1;
    hid_t      under_dapl, vol_id, dacc_tpl;
    hsize_t    my_data_size[DIM], my_offset[DIM], my_count[DIM];
    hsize_t    x_min, x_max, y_min, y_max, z_min, z_max;
    double     t0, t1, t2, t3, t4;
    

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

    opterr = 0;
    strncpy(filename,  "./data/fake_xyz_373G.h5", NAME_MAX);
    strncpy(global_query, "[100:200,300:400,500:600]", NAME_MAX);
    //strncpy(global_query, "[:,:,20000:40000]", NAME_MAX);
    strncpy(dataset_name,"/entry_0/data_0/data_0", NAME_MAX);

    while ((c = getopt (argc, argv, "f:d:q:c:")) != -1) {
        switch (c)
        {
        case 'f':
            strncpy(filename, optarg, NAME_MAX);
            break;
        case 'q':
            strncpy(global_query, optarg, NAME_MAX);
            break;
        case 'd':
            strncpy(dataset_name, optarg, NAME_MAX);
            break;
        case 'c':
            cio = atoi(optarg);
            break;
        default:
            printf("Error option [%s]\n", optarg);
            exit(-1);
        }
    }

    ht_query_parser(global_query, &x_min, &x_max, &y_min, &y_max, &z_min, &z_max);

    MPI_Barrier(MPI_COMM_WORLD);
    t0 = MPI_Wtime();


    /* Open the file */
    hid_t    file_id_2, dataset_id_2, dset_space_id_2, plist_id;

    plist_id = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL);

    file_id_2       = H5Fopen(filename, H5F_ACC_RDONLY, plist_id);
    dataset_id_2    = H5Dopen(file_id_2, dataset_name, H5P_DEFAULT);
    dset_space_id_2 = H5Dget_space(dataset_id_2);

    // Get dataset dimension
    H5Sget_simple_extent_dims(dset_space_id_2, my_data_size, NULL);

    // Change to max dimension of dataset if applies
    if(x_max == -1)
        x_max = my_data_size[0];
    if(y_max == -1)
        y_max = my_data_size[1];
    if(z_max == -1)
        z_max = my_data_size[2];

    MPI_Barrier(MPI_COMM_WORLD);
    t1 = MPI_Wtime();

    if (mpi_rank == 0) {
        printf("Open file: [%s] and dataset [%s], takes [%f]s ", filename, "/entry_0/data_0/data_0", (t1-t0));
        printf("[%llu:%llu,%llu:%llu,%llu:%llu]\n", x_min, x_max, y_min, y_max, z_min, z_max);
    }


    my_offset[0] = x_min;
    my_offset[1] = y_min;
    my_offset[2] = z_min + ((z_max - z_min)/mpi_size*mpi_rank);

    my_count[0] = x_max - x_min;
    my_count[1] = y_max - y_min;
    if(mpi_rank != (mpi_size -1)){
        my_count[2] = (z_max - z_min)/mpi_size;
    }
    else{
        my_count[2] = (z_max - z_min)/mpi_size +  (z_max - z_min)%mpi_size;
    }

    // when z is less than proc_num
    // TODO: this is a hard code for 1024 proc
    //       need correction!
    if(my_count[2] == 0 || my_count[2] == (z_max - z_min)) {
        if(mpi_rank %2 == 1) {
            my_offset[1] += (y_max-y_min)/2;
        }

        my_count[1] = (y_max - y_min)/2;
        my_count[2] = 1;
    }

    hid_t    memspace;
    memspace =  H5Screate_simple(3, my_count, NULL);

    H5Sselect_hyperslab(dset_space_id_2, H5S_SELECT_SET, my_offset, NULL, my_count, NULL);
    if (mpi_rank == 0 || mpi_rank == (mpi_size -1)) {
        printf("\n%d: X Y Z start = [%llu, %llu, %llu]\n", mpi_rank,  my_offset[0], my_offset[1], my_offset[2]);
        printf("%d: X Y Z size  = [%llu, %llu, %llu]\n", mpi_rank, my_count[0], my_count[1], my_count[2]);
    }

    hsize_t         buff_size = my_count[0] * my_count[1] * my_count[2];
    unsigned short *data_buf2 = (unsigned short *)malloc(buff_size * sizeof(unsigned short));
    if(data_buf2 == NULL) {
        printf("Memory allication fails ! \n");
    }

    if(cio == 1) {
        hid_t plist2_id = H5Pcreate(H5P_DATASET_XFER);
        H5Pset_dxpl_mpio(plist2_id, H5FD_MPIO_COLLECTIVE);
        H5Dread(dataset_id_2, H5T_NATIVE_USHORT, memspace, dset_space_id_2, plist2_id, data_buf2);
        H5Pclose(plist2_id);
    } else {
        H5Dread(dataset_id_2, H5T_NATIVE_USHORT, memspace, dset_space_id_2, H5P_DEFAULT, data_buf2);
    }

    MPI_Barrier(MPI_COMM_WORLD);
    t3 = MPI_Wtime();
    if(mpi_rank ==0 ) {
        printf("Rank %llu (%llu, %llu, %llu), Reading takes [%f]s \n", mpi_rank, my_data_size[0], my_data_size[1], my_data_size[2], (t3-t1));
    }

    H5Sclose(dset_space_id_2);
    H5Dclose(dataset_id_2);
    H5Pclose(plist_id);
    H5Fclose(file_id_2);

    MPI_Barrier(MPI_COMM_WORLD);
    t4 = MPI_Wtime();
    if(mpi_rank ==0 ) {
        printf("File and dset closing takes [%f]s,  Overall takes [%f]s \n", t4-t3, t4-t0);
    }

    if(data_buf2 != NULL)
        free(data_buf2);

    H5close();
    MPI_Finalize();
    return 0;
}

