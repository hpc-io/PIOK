#include <stdio.h>
#include <stdlib.h>
#include <hdf5.h>

#include "params.h"
#include "read_h5.h"

// Assuming each time step data refers to 1 day's data
// This version reads a block of 2D vapor array using start [] and count []
//      values that are specific to a desired region
//      from one time step
long int read_vapor_data_nc_1day (float * vapor_data, input_params * ips, \
                                  int rows_from, int rows_to,           \
                                  int cols_from, int cols_to)
{
    int retval;
    int pos, start_pos;
    char dset_name[128];
    int i,j;

    hid_t file_id, space_id, memspace_id, dataset_id, plist_id;

    // Find the block size to read
    int num_cols = cols_to - cols_from;
    int num_rows = rows_to - rows_from;

    // Allocate buffer for partial data for each time step
    float part_vapor_data[1][num_rows][num_cols];

    // printf ("File name: %s \n", ips->ifname[2]);
    // init vapor data to 0 -- Just to make sure there is no garbage
    for ( i = 0; i < num_rows; i++)
    {
        for ( j = 0; j < num_cols; j++)
        {
            pos = i * num_cols + j;
            vapor_data[pos] = 0;
        }
    }

    // dataset path
    sprintf(dset_name, "/%s", DATA_VAR_NAME);

    // Open 1 files; Assuming that ifname[2] has the filename
    plist_id    = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL);

    /* printf("%s:%s\n", ips->ifname[2], dset_name); */

    file_id     = H5Fopen(ips->ifname[2], H5F_ACC_RDONLY, plist_id);
    dataset_id  = H5Dopen(file_id, dset_name, H5P_DEFAULT);
    space_id    = H5Dget_space(dataset_id);

    int iter = 2;
    start_pos = ips->time_step[iter];

    // printf ("start_pos: %d rows_from: %d cols_from: %d num_rows: %d  num_cols: %d \n", start_pos, rows_from, cols_from, num_rows, num_cols);

    hsize_t start[3] = {start_pos, rows_from, cols_from};
    hsize_t count[3] = {1, num_rows, num_cols};

    /* hsize_t data_size[2]; */
    /* data_size[0] = num_rows; */
    /* data_size[0] = num_cols; */

    memspace_id = H5Screate_simple(3, count, NULL); 

    H5Sselect_hyperslab(space_id, H5S_SELECT_SET, start, NULL, count, NULL);

    H5Dread(dataset_id, H5T_IEEE_F32LE, memspace_id, space_id,
                         H5P_DEFAULT, part_vapor_data);


    for (i = 0; i < num_rows; i++) 
    {    
       for (j = 0; j < num_cols; j++) 
       {    
           pos = i * num_cols + j; 
           vapor_data[pos] = part_vapor_data[0][i][j];
           // printf ("%.3f   ", part_vapor_data[i][j]);
       }    
       // printf ("\n");
    }  


    H5Dclose(dataset_id);
    H5Sclose(space_id);
    H5Sclose(memspace_id);
    H5Fclose(file_id);


    return (num_rows * num_cols);
}
// end read_vapor_data_nc_1day ()

