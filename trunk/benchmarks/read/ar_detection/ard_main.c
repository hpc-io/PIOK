/**********************
ard_main.cpp
Author: Suren Byna <SByna@lbl.gov>
Lawrence Berkeley Lab

Description:
This is an implementation to detect Atmospheric Rivers in a given vapor data. Current implementation reads vapor data from a given file. The ar_detect funtion thresholds the vapor data to a low and high values and then labels the data. The labeled data is then passed to SAUF (Scan with Array-based Union Find) funtion to identify connected  components in 2D array. The generated connected component labels are then scanned to verify whether any vapor bands (atmospheric rivers) originating tropical area and ending at any land. When high moisture touches land, there could be a lot of rain fall, which would affect human lives. The vapor thresholds, region of interest, and origin and destination of vapor bands are parameterized to detect atmospheric rivers in any given region of the world map.

**********************/

//#include "ar_detect.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <unistd.h>
#include <math.h>
#include "mpi.h"

#include "params.h"
#include "read_h5.h"
#include <hdf5.h>

#define PROFILE

#define NUM_THR_RANGES 4
#define pi 3.14159265358979323846

#ifdef PROFILE
#include "timer.h"
struct timeval start_time[15];
float elapse[15];
#endif

//parse_args: Parse arguments
static const char *options = "f:i:o:s:";

int parse_args (int argc, char **argv, input_params *ips)
{
	opterr = 0;
	int index;
	int c;

	while ((c = getopt (argc, argv, options)) != -1)
	{
		switch (c)
		{
			case 'f':
				ips->file_type = atoi (optarg);
				break;
			case 'i':
				ips->ifname[0] = (optarg);
				ips->ifname[1] = (optarg);
				ips->ifname[2] = (optarg);
				break;
			case 'o':
				ips->ofname = optarg;
				break;
			case 's':
				ips->dataset = optarg;
				break;
			default:
				return -1;
		}
	}

	for (index = optind; index < argc; index++)
        {
                printf ("Non-option argument %s\n", argv[index]);
                return -1;
        }
	if (ips->file_type == 0)
	{
		printf ("Invalid file_type, options 1 or 2. \n");
		return -1;
	}
	if (ips->ifname == NULL)
	{
		printf ("Empty input file name not allowed. \n");
		return -1;
	}

	if (ips->ofname == NULL)
	{
		printf ("Empty output file name not allowed. \n");
		return -1;
	}

	return 1;
}

void print_usage ()
{
	printf ("\tUsage: ar_detect -f input_type (1 or 2) \n");
	printf ("\t\t-i input_file_name -o output_file_name -s dataset_type (CAM5 or CMIP5\n");
	printf ("\n");
}

void print_input_params (input_params *ips)
{
	printf ("Input Parameters: \n");
	printf ("\tNumber of lats: %d \n", ips->num_lats);
	printf ("\tNumber of lons: %d \n", ips->num_lons);
	if (ips->file_type == 1)
	{
		printf ("\tFile type: One input file\n");
	}
	else if (ips->file_type == 2)
	{
		printf ("\tFile type: Input of input files\n");
	}
	printf ("\tInput file name 1: %s \n", ips->ifname[0]);
	printf ("\tInput file name 2: %s \n", ips->ifname[1]);
	printf ("\tInput file name 3: %s \n", ips->ifname[2]);
	printf ("\tTime step 1: %d\n", ips->time_step[0]);
	printf ("\tTime step 2: %d\n", ips->time_step[1]);
	printf ("\tTime step 3: %d\n", ips->time_step[2]);
	printf ("\tOutput file name: %s \n", ips->ofname);
}

// The main function
int main (int argc, char **argv)
{
	int my_rank;

	input_params *ips = (input_params *) malloc (sizeof (input_params));

	// Parse arguments
	int retval = parse_args (argc, argv, ips);
	if (retval < 0)
	{
		print_usage ();
		free (ips);
		exit (1);
	}
	
	// Init MPI
	MPI_Init (&argc, &argv);
	MPI_Comm_rank (MPI_COMM_WORLD, &my_rank);

	timer_on (1);

	// std::vector< std::string > filenames(3);
	// std::vector< int > time_steps (3);
	/* vector< string > filenames(3); */
	/* vector< int > time_steps (3); */

	/* filenames[0] = ips->ifname[0]; */
	/* filenames[1] = ips->ifname[1]; */
	/* filenames[2] = ips->ifname[2]; */
	
	/* if (my_rank == 0) */
	/* { */
	/* 	time_steps[0] = my_rank; */
	/* 	time_steps[1] = my_rank; */
	/* 	time_steps[2] = my_rank; */
	/* } */
	/* else if (my_rank == 1) */
	/* { */
	/* 	time_steps[0] = my_rank; */
	/* 	time_steps[1] = my_rank-1; */
	/* 	time_steps[2] = my_rank-1; */
	/* } */
	/* else if (my_rank >= 2) */
	/* { */
	/* 	time_steps[0] = my_rank; */
	/* 	time_steps[1] = my_rank-1; */
	/* 	time_steps[2] = my_rank-2; */
	/* } */

	// The driver function
	//ar_detect_nc (filenames, time_steps, ips->ofname, ips->dataset);
        int cols_from, cols_to, rows_from, rows_to;
        cols_from = 0;
        cols_to   = 16;
        rows_from = 0;
        rows_to   = 16;

        float* data;
        data = (float*)malloc(sizeof(float)*(cols_to-cols_from)*(rows_to-rows_from));

        read_vapor_data_nc_1day(data, ips, rows_from, rows_to, cols_from, cols_to);

	timer_off (1);

	free (ips);

        H5close();

	MPI_Finalize ();
	return 0;
}
