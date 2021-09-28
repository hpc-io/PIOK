/****** Copyright Notice ***
 *
 * PIOK - Parallel I/O Kernels - VPIC-IO, VORPAL-IO, and GCRM-IO, Copyright
 * (c) 2015, The Regents of the University of California, through Lawrence
 * Berkeley National Laboratory (subject to receipt of any required
 * approvals from the U.S. Dept. of Energy).  All rights reserved.
 *
 * If you have questions about your rights to use or distribute this
 * software, please contact Berkeley Lab's Innovation & Partnerships Office
 * at  IPO@lbl.gov.
 *
 * NOTICE.  This Software was developed under funding from the U.S.
 * Department of Energy and the U.S. Government consequently retains
 * certain rights. As such, the U.S. Government has been granted for itself
 * and others acting on its behalf a paid-up, nonexclusive, irrevocable,
 * worldwide license in the Software to reproduce, distribute copies to the
 * public, prepare derivative works, and perform publicly and display
 * publicly, and to permit other to do so.
 *
 ****************************/

/**
 *
 * Email questions to SByna@lbl.gov
 * Scientific Data Management Research Group
 * Lawrence Berkeley National Laboratory
 *
*/


// Description: This is a simple benchmark based on VPIC's I/O interface
//		Each process writes a specified number of particles into 
//		a hdf5 output file using H5Dwrite_multi (), where HDF5 writes 
//		multiple datasets with one call
//		create_and_write_synthetic_h5_data () function writes datasets 
//			using seperate H5Dwrite () calls
//		create_and_write_synthetic_h5_md_data () function writes datasets 
//			using H5Dwrite_multi () call
// Author:	Suren Byna <SByna@lbl.gov>
//		Lawrence Berkeley National Laboratory, Berkeley, CA
// Created:	in 2011
// Modified:	01/06/2014 --> Removed all H5Part calls and using HDF5 calls
// Modified:	03/12/2014 --> Added create_and_write_synthetic_h5_md_data () function
//			       for writing multiple HDF5 datasets using one 
//			       H5Dwrite_multi () call
// 

#include <assert.h>
#include "hdf5.h"
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>

#define MULTI_DSET 1

// A simple timer based on gettimeofday
#include "./timer.h"
struct timeval start_time[8];
float elapse[8];

// HDF5 specific declerations
herr_t ierr;
hid_t file_id, dset_id;
hid_t filespace, memspace;
hid_t plist_id;

// Variables and dimensions
long numparticles = 8388608;	// 8  meg particles per process
long long total_particles, offset;
#define NUM_DSETS 8
#define FAIL -1

float *x, *y, *z;
float *px, *py, *pz;
int *id1, *id2;
int x_dim = 64;
int y_dim = 64; 
int z_dim = 64;

// Uniform random number
inline double uniform_random_number() 
{
	return (((double)rand())/((double)(RAND_MAX)));
}

// Initialize particle data
void init_particles ()
{
	int i;
	for (i=0; i<numparticles; i++) 
	{
		id1[i] = i;
		id2[i] = i*2;
		x[i] = uniform_random_number()*x_dim;
		y[i] = uniform_random_number()*y_dim;
		z[i] = ((double)id1[i]/numparticles)*z_dim;    
		px[i] = uniform_random_number()*x_dim;
		py[i] = uniform_random_number()*y_dim;
		pz[i] = ((double)id2[i]/numparticles)*z_dim;    
	}
}

// Create HDF5 file and write data using different H5Dwrite calls
void create_and_write_synthetic_h5_data(int rank)
{
	// Note: printf statements are inserted basically 
	// to check the progress. Other than that they can be removed
	dset_id = H5Dcreate(file_id, "x", H5T_NATIVE_FLOAT, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        ierr = H5Dwrite(dset_id, H5T_NATIVE_FLOAT, memspace, filespace, plist_id, x);
        H5Dclose(dset_id);
	if (rank == 0) printf ("Written variable 1 \n");

	dset_id = H5Dcreate(file_id, "y", H5T_NATIVE_FLOAT, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        ierr = H5Dwrite(dset_id, H5T_NATIVE_FLOAT, memspace, filespace, plist_id, y);
        H5Dclose(dset_id);
	if (rank == 0) printf ("Written variable 2 \n");

	dset_id = H5Dcreate(file_id, "z", H5T_NATIVE_FLOAT, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        ierr = H5Dwrite(dset_id, H5T_NATIVE_FLOAT, memspace, filespace, plist_id, z);
        H5Dclose(dset_id);
	if (rank == 0) printf ("Written variable 3 \n");

	dset_id = H5Dcreate(file_id, "id1", H5T_NATIVE_INT, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        ierr = H5Dwrite(dset_id, H5T_NATIVE_INT, memspace, filespace, plist_id, id1);
        H5Dclose(dset_id);
	if (rank == 0) printf ("Written variable 4 \n");

	dset_id = H5Dcreate(file_id, "id2", H5T_NATIVE_INT, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        ierr = H5Dwrite(dset_id, H5T_NATIVE_INT, memspace, filespace, plist_id, id2);
        H5Dclose(dset_id);
	if (rank == 0) printf ("Written variable 5 \n");

	dset_id = H5Dcreate(file_id, "px", H5T_NATIVE_FLOAT, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        ierr = H5Dwrite(dset_id, H5T_NATIVE_FLOAT, memspace, filespace, plist_id, px);
        H5Dclose(dset_id);
	if (rank == 0) printf ("Written variable 6 \n");

	dset_id = H5Dcreate(file_id, "py", H5T_NATIVE_FLOAT, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        ierr = H5Dwrite(dset_id, H5T_NATIVE_FLOAT, memspace, filespace, plist_id, py);
        H5Dclose(dset_id);
	if (rank == 0) printf ("Written variable 7 \n");

	dset_id = H5Dcreate(file_id, "pz", H5T_NATIVE_FLOAT, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        ierr = H5Dwrite(dset_id, H5T_NATIVE_FLOAT, memspace, filespace, plist_id, pz);
        H5Dclose(dset_id);
	if (rank == 0) printf ("Written variable 8 \n");
}

// Create and write datasets using one multi-dataset write
void create_and_write_synthetic_h5_md_data(int rank)
{
	H5D_rw_multi_t dset_info[NUM_DSETS];
	hid_t dset_id_x, dset_id_y, dset_id_z; 
	hid_t dset_id_id1, dset_id_id2;
	hid_t dset_id_px, dset_id_py, dset_id_pz;
	
	//init dset_info[]
	int i;
	for (i=0; i<NUM_DSETS; i++) 
	{
		memset(&dset_info[i], 0, sizeof(H5D_rw_multi_t));
	}
	
	//create datasets collectively
	dset_id_x = H5Dcreate(file_id, "x", H5T_NATIVE_FLOAT, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	assert(dset_id_x != FAIL);

	dset_id_y = H5Dcreate(file_id, "y", H5T_NATIVE_FLOAT, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	assert(dset_id_y != FAIL);

	dset_id_z = H5Dcreate(file_id, "z", H5T_NATIVE_FLOAT, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	assert(dset_id_z != FAIL);

	dset_id_id1 = H5Dcreate(file_id, "id1", H5T_NATIVE_INT, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	assert(dset_id_id1 != FAIL);

	dset_id_id2 = H5Dcreate(file_id, "id2", H5T_NATIVE_INT, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	assert(dset_id_id2 != FAIL);

	dset_id_px = H5Dcreate(file_id, "px", H5T_NATIVE_FLOAT, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	assert(dset_id_px != FAIL);

	dset_id_py = H5Dcreate(file_id, "py", H5T_NATIVE_FLOAT, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	assert(dset_id_py != FAIL);

	dset_id_pz = H5Dcreate(file_id, "pz", H5T_NATIVE_FLOAT, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	assert(dset_id_pz != FAIL);


	//set up dset_info structure
	dset_info[0].dset_id = dset_id_x;
	dset_info[0].mem_type_id = H5T_NATIVE_FLOAT;
	dset_info[0].mem_space_id = memspace;
	dset_info[0].dset_space_id = filespace;
	dset_info[0].u.wbuf = x;
	
	dset_info[1].dset_id = dset_id_y;
	dset_info[1].mem_type_id = H5T_NATIVE_FLOAT;
	dset_info[1].mem_space_id = memspace;
	dset_info[1].dset_space_id = filespace;
	dset_info[1].u.wbuf = y;
	
	dset_info[2].dset_id = dset_id_z;
	dset_info[2].mem_type_id = H5T_NATIVE_FLOAT;
	dset_info[2].mem_space_id = memspace;
	dset_info[2].dset_space_id = filespace;
	dset_info[2].u.wbuf = z;
	
	dset_info[3].dset_id = dset_id_id1;
	dset_info[3].mem_type_id = H5T_NATIVE_INT;
	dset_info[3].mem_space_id = memspace;
	dset_info[3].dset_space_id = filespace;
	dset_info[3].u.wbuf = id1;
	
	dset_info[4].dset_id = dset_id_id2;
	dset_info[4].mem_type_id = H5T_NATIVE_INT;
	dset_info[4].mem_space_id = memspace;
	dset_info[4].dset_space_id = filespace;
	dset_info[4].u.wbuf = id2;
	
	dset_info[5].dset_id = dset_id_px;
	dset_info[5].mem_type_id = H5T_NATIVE_FLOAT;
	dset_info[5].mem_space_id = memspace;
	dset_info[5].dset_space_id = filespace;
	dset_info[5].u.wbuf = px;
	
	dset_info[6].dset_id = dset_id_py;
	dset_info[6].mem_type_id = H5T_NATIVE_FLOAT;
	dset_info[6].mem_space_id = memspace;
	dset_info[6].dset_space_id = filespace;
	dset_info[6].u.wbuf = py;
	
	dset_info[7].dset_id = dset_id_pz;
	dset_info[7].mem_type_id = H5T_NATIVE_FLOAT;
	dset_info[7].mem_space_id = memspace;
	dset_info[7].dset_space_id = filespace;
	dset_info[7].u.wbuf = pz;
	
	// write multiple dsets
	ierr = H5Dwrite_multi (file_id, plist_id, NUM_DSETS, dset_info);
	if (rank == 0)	{
		assert(ierr != FAIL);
		printf ("Finished writing multiple datasets... \n");
	}
	
	//debugging purpose: write datasets differently using dset_info
	/* 
	for (i = 0; i < NUM_DSETS; i++)
	{
		H5Dwrite(dset_info[i].dset_id, dset_info[i].mem_type_id, dset_info[i].mem_space_id, \
			 dset_info[i].dset_space_id, plist_id, dset_info[i].u.wbuf);
		if (rank == 0)
			printf ("Finished writing variable %d \n", i);
	}
	*/

	//Close datasets
	if (rank == 0) {
		printf ("Before closing all datasets... \n");
	}

	for (i = 0; i < NUM_DSETS; i++)
	{
		if (dset_info[i].dset_id > 0)
		{
			ierr = H5Dclose (dset_info[i].dset_id);
			if (rank == 0)	{
				assert(ierr != FAIL);
			}
		}
	}
	if (rank == 0) {
		printf ("Finished closing all datasets... \n");
	}
}

int main (int argc, char* argv[]) 
{
	char *file_name = argv[1];
	
	MPI_Init(&argc,&argv);
	int my_rank, num_procs;
	MPI_Comm_rank (MPI_COMM_WORLD, &my_rank);
	MPI_Comm_size (MPI_COMM_WORLD, &num_procs);

	MPI_Comm comm  = MPI_COMM_WORLD;
        MPI_Info info  = MPI_INFO_NULL;

	if (argc == 3)
	{
		numparticles = (atoi (argv[2]))*1024*1024;
	}
	else
	{
		numparticles = 8*1024*1024;
	}

	if (my_rank == 0) {printf ("Number of paritcles: %ld \n", numparticles);}

	x=(float*)malloc(numparticles*sizeof(double));
	y=(float*)malloc(numparticles*sizeof(double));
	z=(float*)malloc(numparticles*sizeof(double));

	px=(float*)malloc(numparticles*sizeof(double));
	py=(float*)malloc(numparticles*sizeof(double));
	pz=(float*)malloc(numparticles*sizeof(double));

	id1=(int*)malloc(numparticles*sizeof(int));
	id2=(int*)malloc(numparticles*sizeof(int));

	init_particles ();

	if (my_rank == 0)
	{
		printf ("Finished initializeing particles \n");
	}

	// h5part_int64_t alignf = 8*1024*1024;

	MPI_Barrier (MPI_COMM_WORLD);
	timer_on (0);

	MPI_Allreduce(&numparticles, &total_particles, 1, MPI_LONG_LONG, MPI_SUM, comm);
        MPI_Scan(&numparticles, &offset, 1, MPI_LONG_LONG, MPI_SUM, comm);	
	offset -= numparticles;

	plist_id = H5Pcreate(H5P_FILE_ACCESS);
        H5Pset_fapl_mpio(plist_id, comm, info);
	
	H5Pset_libver_bounds(plist_id, H5F_LIBVER_LATEST, H5F_LIBVER_LATEST);
	H5Pset_link_phase_change(plist_id, 12, 10);

	//file = H5PartOpenFileParallel (file_name, H5PART_WRITE | H5PART_VFD_MPIPOSIX | H5PART_FS_LUSTRE, MPI_COMM_WORLD);
	file_id = H5Fcreate(file_name , H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
        H5Pclose(plist_id);
	
	if (my_rank == 0)
	{
		printf ("Opened HDF5 file... \n");
	}
	// Throttle and see performance
	// H5PartSetThrottle (file, 10);

	// H5PartWriteFileAttribString(file, "Origin", "Tested by Suren");

	filespace = H5Screate_simple(1, (hsize_t *) &total_particles, NULL);
        memspace =  H5Screate_simple(1, (hsize_t *) &numparticles, NULL);

        plist_id = H5Pcreate(H5P_DATASET_XFER);
        // H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
        H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
        // H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_INDEPENDENT);
        H5Sselect_hyperslab(filespace, H5S_SELECT_SET, (hsize_t *) &offset, NULL, (hsize_t *) &numparticles, NULL);

	MPI_Barrier (MPI_COMM_WORLD);
	timer_on (1);

	if (my_rank == 0) printf ("Before writing particles \n");
	if (MULTI_DSET)
	{
		create_and_write_synthetic_h5_md_data(my_rank);
	}
	else 
	{
		create_and_write_synthetic_h5_data(my_rank);
	}

	MPI_Barrier (MPI_COMM_WORLD);
	timer_off (1);
	if (my_rank == 0) printf ("After writing particles \n");

	timer_on (2);
	H5Sclose(memspace);
	if (my_rank == 0) printf ("Closed memspace ...\n");
	timer_off (2);
	
	timer_on (3);
        H5Sclose(filespace);
	if (my_rank == 0) printf ("Closed filespace ...\n");
	timer_off (3);

	timer_on (4);
        H5Pclose(plist_id);
	if (my_rank == 0) printf ("Closed plist_id ...\n");
	timer_off (4);

	timer_on (5);
        H5Fclose(file_id);
	if (my_rank == 0) printf ("After closing HDF5 file \n");
	timer_off (5);

	timer_on (6);
	free(x); free(y); free(z);
	free(px); free(py); free(pz);
	free(id1);
	free(id2);
	timer_off (6);

	MPI_Barrier (MPI_COMM_WORLD);

	timer_off (0);

	if (my_rank == 0)
	{
		printf ("\nTiming results\n");
		timer_msg (1, "just writing data");
		timer_msg (0, "opening, writing, closing file");
		timer_msg (2, "closing memspace");
		timer_msg (3, "closing file space");
		timer_msg (4, "closing property list");
		timer_msg (5, "closing file");
		timer_msg (6, "freeing memory");
		printf ("\n");
	}

	MPI_Finalize();

	return 0;
}
