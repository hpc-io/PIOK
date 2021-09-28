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

/* 
 * File:        gcrmio.c
 * Author:      Mark Howison
 * Created:     2009-05-04
 * Description: Benchmark application that simulates I/O for David Randall's
 *              INCITE19 GCRM code.
 */

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include <H5Part.h>

#define ONE_MEGABYTE 1048576

#define ERROR(msg,my_rank) do{\
    fprintf (stderr, "rank %d: %s (%d): %s!\n",\
            my_rank, __FILE__, __LINE__, msg);\
    fflush (stderr);\
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);\
}while(0);


/* Parameters *****************************************************************/
struct Params {
    /* MPI */
    int         rank;               // MPI rank of this proc
    int         nprocs;             // MPI # of procs

    /* H5Part/HDF5 */
    int         alignment;          // HDF5 alignment in bytes
    int         chunking;           // enable HDF5 chunking
    int         throttle;           // H5Part throtle factor
    int         lustre;             // enable optimizations for lustre
    int         mpiposix;           // enable HDF5 MPI-POSIX VFD
    int         cb;                 // collective buffering

    /* file/grid parameters */
    char*       output;             // path to the test file
    int         t;                  // timesteps
    int         r;                  // grid resolution
    int         s;                  // subdomain resolution
    int         km;                 // number of vertical levels

    /* runtime options */
    int         i;                  // # of iterations
    int         prs;                // write pressure var
    int         edge;               // write edge var
    int         verbosity;          // verbosity level

    /* internal calculations */
    int         nsd;                // number of subdomains
    int         im,jm;              // max indices for rectangular grid coords
    int         nm;                 // subdomains per proc
    long        ngrid;              // size of grid

    /* H5Part state */
    h5part_int64_t ncells;          // total number of cells
    h5part_int64_t ndata;           // total data per task
    h5part_int64_t chunk[3];        // chunk dimensions
    h5part_int64_t layout[6];       // H5Block layout

    /* collective buffering */
    int         isWriter;
    int         nBuffered;
    int         myWriter;
    MPI_Comm    commWriters;
};
typedef struct Params Params;


/******************************************************************************/
/* Argument helper functions                                                  */
/******************************************************************************/

void print_usage (int rank)
{
    if (rank == 0) {
        printf ("Usage: (-r and -s are required!)\n\n");
        printf (" file/grid parameters:\n");
        printf (" -r n\tgrid resolution\n");
        printf (" -s n\tsubdomain resolution\n");
        printf (" -k n\tvertical levels (default 25)\n");
        printf (" -t n\ttimesteps (default 1)\n");
        printf (" -o 's'\toutput path (default 'output')\n\n");
        printf (" runtime options:\n");
        printf (" -i n \titerate n times (default 1)\n");
        printf (" -v n \tverbosity level (default 1)\n");
        printf (" -prs \twrite pressure variable (default off)\n");
        printf (" -edge\twrite edge-centered variable (default off)\n\n");
        printf (" H5Part/HDF5 options:\n");
        printf (" -align n\talign to n byte boundaries (default 1048576)\n");
        printf (" -chunk  \tenable chunking (default off)\n");
        printf (" -lustre \tenable optimizations for lustre (default off)\n");
        printf (" -mpiposix\tenable HDF5 MPI-POSIX VFD (default off)\n");
        printf (" -throttle n\tthrottle the number of concurrent writes (default 0 = off)\n");
        printf (" -cb n\tuse collective buffering with n aggregators (requires mpiposix)\n");
        printf ("\n");
    }
    MPI_Finalize();
    exit (EXIT_SUCCESS);
}

void parse_args (int argc, char** argv, Params* params)
{
    int i;
    int check = 0;

    // default values
    params->t           = 1;
    params->alignment   = ONE_MEGABYTE;
    params->chunking    = 0;
    params->throttle    = 0;
    params->lustre      = 0;
    params->mpiposix    = 0;
    params->cb          = 0;
    params->km          = 25;
    params->output      = "output";
    params->i           = 1;
    params->prs         = 0;
    params->edge        = 0;
    params->verbosity   = 1;
   
    i = 1;
    while (i < argc)
    {
        if (strcmp(argv[i],"-r") == 0)
        {
            i++;
            params->r = atoi (argv[i]);
            check++;
        }
        else if (strcmp(argv[i],"-s") == 0)
        {
            i++;
            params->s = atoi (argv[i]);
            check++;
        }
       // # of vertical levels
        else if (strcmp(argv[i],"-k") == 0)
        {
            i++;
            params->km = atoi (argv[i]);
        }
        else if (strcmp(argv[i],"-t") == 0)
        {
            i++;
            params->t = atoi (argv[i]);
        }
        else if (strcmp(argv[i],"-i") == 0)
        {
            i++;
            params->i = atoi (argv[i]);
        }
        else if (strcmp(argv[i],"-v") == 0)
        {
            i++;
            params->verbosity = atoi (argv[i]);
        }
        else if (strcmp(argv[i],"-o") == 0)
        {
            i++;
            params->output = (char*) malloc (strlen (argv[i]) + 1);
            strcpy (params->output, argv[i]);
        }
        else if (strcmp(argv[i],"-prs") == 0)
        {
            params->prs = 1;
        }
        else if (strcmp(argv[i],"-edge") == 0)
        {
            params->edge = 1;
        }
        else if (strcmp(argv[i],"-align") == 0)
        {
            i++;
            params->alignment = atoi (argv[i]);
        }
        else if (strcmp(argv[i],"-chunk") == 0)
        {
            params->chunking = 1;
        }
        else if (strcmp(argv[i],"-lustre") == 0)
        {
            params->lustre = 1;
        }
        else if (strcmp(argv[i],"-mpiposix") == 0)
        {
            params->mpiposix = 1;
        }
        else if (strcmp(argv[i],"-throttle") == 0)
        {
            i++;
            params->throttle = atoi (argv[i]);
        }
        else if (strcmp(argv[i],"-cb") == 0)
        {
            i++;
            params->cb = atoi (argv[i]);
        }
        // print usage
        else if (strcmp(argv[i],"--help") == 0)
        {
            print_usage (params->rank);
        }
        else
        {
            if (params->rank == 0) {
                fprintf (stderr, "%s: unrecognized argument %s \n",
                        argv[0], argv[i]);
            }
            print_usage (params->rank);
        }
        i++;
    }

    if (check != 2) print_usage (params->rank);

    if (params->rank == 0) {
        printf ("Parameters:\n");
        printf ("\tgrid resolution = %d\n", params->r);
        printf ("\tsubdomain resolution = %d\n", params->s);
        printf ("\ttimesteps = %d\n", params->t);
        printf ("\tvertical levels = %d\n", params->km);
        printf ("\toutput path = '%s'\n", params->output);
        printf ("\tverbosity level = %d\n", params->verbosity);
        printf ("\tprs = %d\n", params->prs);
        printf ("\tedge = %d\n", params->edge);
        printf ("\talignment = %d\n", params->alignment);
        printf ("\tchunking = %d\n", params->chunking);
        printf ("\tlustre = %d\n", params->lustre);
        printf ("\tmpiposix = %d\n", params->mpiposix);
        printf ("\tthrottle = %d\n", params->throttle);
        printf ("\tcb = %d\n", params->cb);
    }

    if (params->cb) {
        if (! params->mpiposix) {
            ERROR("Collective buffering only works with MPI-POSIX",
                    params->rank);
        }
        if (params->throttle) {
            ERROR("Collective buffering is incompatible with throttling",
                    params->rank);
        }
        if (params->cb > params->nprocs) {
            ERROR("Fewer MPI tasks than CB aggregators", params->rank);
        }
    }

    if (params->r < 5 || params->r > 16) {
        ERROR("Grid resolution must be in the range [5..16]",
                params->rank);
    }

    if (params->s < 0 || params->s > 7) {
        ERROR("Subdomain resolution must be in the range [0..7]",
                params->rank);
    }

    if (params->r <= params->s) {
        ERROR("Grid resolution must be larger than subdomain resolution",
                params->rank);
    }

    // compute # of subdomains
    params->nsd = 10;
    for (i=0;i<params->s;i++) params->nsd *= 4;

    if (params->rank == 0) {
        printf ("Computed values:\n");
        printf ("\tsubdomains = %d\n", params->nsd);
    }

    if (params->nsd % params->nprocs != 0) {
        ERROR("Number of procs must divide number of subdomains",
                params->rank);
    }

    // compute 2D extents
    params->im = 1;
    for (i=0;i<(params->r - params->s);i++) params->im *= 2;
    params->jm = params->im;

    // compute size of grid
    params->ngrid = 10;
    for (i=0;i<params->r;i++) params->ngrid *= 4;
    // add two for the north and south poles
    params->ngrid += 2;

    // compute subdomains per proc
    params->nm = params->nsd / params->nprocs;

    if (params->rank == 0) {
        double chunk, size;
        printf ("\tgrid size = %ld\n", params->ngrid);
        printf ("\ti = %d\n", params->im);
        printf ("\tj = %d\n", params->jm);
        printf ("\tn = %d\n", params->nm);
        chunk  = params->im;
        chunk *= params->jm;
        chunk *= params->km + 1;
        chunk *= sizeof(float);
        chunk /= 1024;
        printf ("\tprs chunk size = %g KB\n", chunk);
        size  = params->ngrid;
        size *= params->km+1;
        size *= sizeof(float);
        size /= ONE_MEGABYTE;
        printf ("\tprs timestep size = %g MB\n", size);
        size *= params->t;
        printf ("\tprs file size = %g MB\n", size);
        chunk  = params->im;
        chunk *= params->jm;
        chunk *= params->km;
        chunk *= 6.0;
        chunk *= sizeof(float);
        chunk /= 1024;
        printf ("\tedge chunk size = %g KB\n", chunk);
        size  = params->ngrid;
        size *= params->km;
        size *= 6.0;
        size *= sizeof(float);
        size /= ONE_MEGABYTE;
        printf ("\tedge timestep size = %g MB\n", size);
        size *= params->t;
        printf ("\tedge file size = %g MB\n", size);
    }
}


/******************************************************************************/
/* Write functions                                                            */
/******************************************************************************/

void init_cb1(Params *params)
{
    // the bottom #cb ranks are the writers
    if (params->rank / params->cb == 0) {
        params->isWriter = 1;
    } else {
        params->isWriter = 0;
    }

    params->myWriter = params->rank % params->cb;
    params->nBuffered = (params->nprocs - params->myWriter - 1) / params->cb;

    int status;
    int range[1][3];
    MPI_Group MPI_GROUP_WORLD;
    MPI_Group group;

    MPI_Comm_group(MPI_COMM_WORLD, &MPI_GROUP_WORLD);

    range[0][0] = 0; // start
    range[0][1] = params->cb - 1; // end
    range[0][2] = 1; // sride

    status = MPI_Group_range_incl(MPI_GROUP_WORLD, 1, range, &group);
    if (status != MPI_SUCCESS) ERROR("MPI_Group_incl failed", params->rank);

    status = MPI_Comm_create(MPI_COMM_WORLD, group, &(params->commWriters));
    if (status != MPI_SUCCESS) ERROR("MPI_Comm_create failed", params->rank);
}

void init_cb2(Params *params)
{
    params->isWriter = 0;

    // writers are spread out across the ranks
    int cbwidth = params->nprocs / params->cb;
    int remainder = params->nprocs - (cbwidth * params->cb);
    int cutoff = remainder * (cbwidth+1);
    if (params->rank < cutoff) {
        if ((params->rank % (cbwidth+1)) == 0) params->isWriter = 1;
        params->myWriter = (cbwidth+1) * (params->rank / (cbwidth+1));
        params->nBuffered = cbwidth + 1;
    } else {
        int tmprank = params->rank - cutoff;
        if (tmprank % cbwidth == 0) params->isWriter = 1;
        params->myWriter = cutoff + cbwidth * (tmprank / cbwidth);
        params->nBuffered = cbwidth;
    }

    int status;
    int range[2][3];
    MPI_Group MPI_GROUP_WORLD;
    MPI_Group group;

    MPI_Comm_group(MPI_COMM_WORLD, &MPI_GROUP_WORLD);

    if (cutoff > 0)
    {
        range[0][0] = 0; // start
        range[0][1] = cutoff - 1; // end
        range[0][2] = cbwidth + 1; // sride

        range[1][0] = cutoff; // start
        range[1][1] = params->nprocs - 1; // end
        range[1][2] = cbwidth; // sride

        status = MPI_Group_range_incl(MPI_GROUP_WORLD, 2, range, &group);
        if (status != MPI_SUCCESS) ERROR("MPI_Group_incl failed", params->rank);
    }
    else
    {
        range[0][0] = 0; // start
        range[0][1] = params->nprocs - 1; // end
        range[0][2] = cbwidth; // sride

        status = MPI_Group_range_incl(MPI_GROUP_WORLD, 1, range, &group);
        if (status != MPI_SUCCESS) ERROR("MPI_Group_incl failed", params->rank);
    }

    status = MPI_Comm_create(MPI_COMM_WORLD, group, &(params->commWriters));
    if (status != MPI_SUCCESS) ERROR("MPI_Comm_create failed", params->rank);
}

void init_cb3(Params *params)
{
    params->isWriter = 0;

    // writers are spread out across the ranks
    int cbwidth = (int)ceil((float)params->nprocs / (float)params->cb);
    if ((cbwidth * params->cb - params->nprocs) >= cbwidth) {
        ERROR("Number of cb nodes is incompatible", params->rank);
    }
    if ((params->rank % cbwidth) == 0) params->isWriter = 1;
    params->myWriter = cbwidth * (params->rank / cbwidth);
    params->nBuffered = cbwidth;
    
    /* debug */
    if (params->isWriter) {
        printf("%d: buffering %d\n", params->rank, params->nBuffered);
    }

    int status;
    int range[2][3];
    MPI_Group MPI_GROUP_WORLD;
    MPI_Group group;

    MPI_Comm_group(MPI_COMM_WORLD, &MPI_GROUP_WORLD);

    {
        range[0][0] = 0; // start
        range[0][1] = params->nprocs - 1; // end
        range[0][2] = cbwidth; // sride

        status = MPI_Group_range_incl(MPI_GROUP_WORLD, 1, range, &group);
        if (status != MPI_SUCCESS) ERROR("MPI_Group_incl failed", params->rank);
    }

    status = MPI_Comm_create(MPI_COMM_WORLD, group, &(params->commWriters));
    if (status != MPI_SUCCESS) ERROR("MPI_Comm_create failed", params->rank);
}

void write_h5part_all (
        H5PartFile *file,
        float* data,
        Params* params)
{
    int i;
    h5part_int64_t status;

    if (params->chunking) {
        status = H5BlockDefine3DChunkDims (file,
                params->chunk[0], params->chunk[1], params->chunk[2]);
        if (status != H5PART_SUCCESS) {
            ERROR("H5BlockDefine3DChunk failed", params->rank);
        }
    }

    status = H5BlockDefine3DFieldLayout (file,
            params->layout[0], params->layout[1],
            params->layout[2], params->layout[3],
            params->layout[4], params->layout[5]);
    if (status != H5PART_SUCCESS) {
        ERROR("H5BlockDefine3DFieldLayout failed", params->rank);
    }

    for (i=0; i<params->t; i++)
    {
        status = H5PartSetStep (file, i);
        if (status != H5PART_SUCCESS) {
            ERROR("H5PartSetStep failed", params->rank);
        }

        status = H5Block3dWriteScalarFieldFloat32 (file, "data", data);
        if (status != H5PART_SUCCESS) {
            ERROR("H5Block3dWriteScalarFieldFloat32 failed",
                    params->rank);
        }
    }
}

void write_h5part_cb1 (
        H5PartFile *file,
        float* data,
        Params* params)
{
    int i;
    double start_time, comm_time;
    h5part_int64_t status;
    int range[1][3];
    MPI_Group MPI_GROUP_WORLD;
    MPI_Group group;
    MPI_Comm comm;
    float *buffer;

    start_time = MPI_Wtime();
    comm_time = 0;

    // create new communicator for aggregation
    range[0][0] = params->myWriter; // start
    range[0][1] = params->myWriter + params->cb*params->nBuffered; // end
    range[0][2] = params->cb; // sride

    MPI_Comm_group(MPI_COMM_WORLD, &MPI_GROUP_WORLD);

    status = MPI_Group_range_incl(MPI_GROUP_WORLD, 1, range, &group);
    if (status != MPI_SUCCESS) ERROR("MPI_Group_range_incl failed", params->rank);

    status = MPI_Comm_create(MPI_COMM_WORLD, group, &comm);
    if (status != MPI_SUCCESS) ERROR("MPI_Comm_create failed", params->rank);

    if (params->isWriter) {
        buffer = (float *) malloc((params->nBuffered+1) * params->ndata * sizeof(float));
        if (!buffer) {
            ERROR("Could not allocate buffer on aggregator node",
                    params->rank);
        }
    }

    status = MPI_Barrier(MPI_COMM_WORLD);

    comm_time += MPI_Wtime() - start_time;
    
    if (params->isWriter)
    {
        status = H5PartSetNumParticles (file,
                (params->nBuffered+1) * params->ndata);
        if (status != H5PART_SUCCESS) {
            ERROR("H5PartSetNumParticles failed", params->rank);
        }

        if (params->chunking) {
            status = H5PartSetChunkSize (file, params->ndata);
            if (status != H5PART_SUCCESS) {
                ERROR("H5PartSetChunkSize failed", params->rank);
            }
        }
    }

    for (i=0; i<params->t; i++)
    {
        // communication phase
        start_time = MPI_Wtime();

        status = MPI_Gather(
                data, (int)params->ndata, MPI_FLOAT,
                buffer, (int)params->ndata, MPI_FLOAT,
                0, comm);
        if (status != MPI_SUCCESS) ERROR("MPI_Gather failed", params->rank);

        comm_time += MPI_Wtime() - start_time;

        if (params->isWriter)
        {
            status = H5PartSetStep (file, i);
            if (status != H5PART_SUCCESS) {
                ERROR("H5PartSetStep failed", params->rank);
            }

            int j;
            for (j=0; j<(params->nBuffered+1); j++)
            {
                h5part_int64_t start = (params->rank + j*params->cb)*params->ndata;
                status = H5PartSetView (file, start, start + params->ndata - 1);
                if (status != H5PART_SUCCESS) {
                    ERROR("H5PartSetView failed", params->rank);
                }

                status = H5PartWriteDataFloat32 (file, "data", buffer);
                if (status != H5PART_SUCCESS) {
                    ERROR("H5PartWriteDataFloat32 failed",
                            params->rank);
                }
            }
        } 
    }

    if (params->isWriter) free(buffer);

    MPI_Group_free(&MPI_GROUP_WORLD);
    MPI_Group_free(&group);
    MPI_Comm_free(&comm);

    if (params->rank == 0) {
        printf ("(communication time = %g s)\n", comm_time);
    }
}

void write_h5part_cb2 (
        H5PartFile *file,
        float* data,
        Params* params)
{
    int i;
    double start_time, comm_time;
    h5part_int64_t status;
    int range[1][3];
    MPI_Group MPI_GROUP_WORLD;
    MPI_Group group;
    MPI_Comm comm;
    float *buffer;

    start_time = MPI_Wtime();
    comm_time = 0;

    int cbwidth = params->nBuffered;
    /*if ((params->myWriter + cbwidth) > params->nprocs)
                cbwidth = params->nprocs - params->rank;*/
    //printf("%d: writer %d cbwidth %d\n", params->rank, params->myWriter, cbwidth);

    // create new communicator for aggregation
    range[0][0] = params->myWriter; // start
    range[0][1] = params->myWriter + cbwidth - 1; // end
    range[0][2] = 1; // sride

    MPI_Comm_group(MPI_COMM_WORLD, &MPI_GROUP_WORLD);

    status = MPI_Group_range_incl(MPI_GROUP_WORLD, 1, range, &group);
    if (status != MPI_SUCCESS) ERROR("MPI_Group_range_incl failed", params->rank);

    status = MPI_Comm_create(MPI_COMM_WORLD, group, &comm);
    if (status != MPI_SUCCESS) ERROR("MPI_Comm_create failed", params->rank);

    if (params->isWriter) {
        buffer = (float *) malloc(params->nBuffered * params->ndata * sizeof(float));
        if (!buffer) {
            ERROR("Could not allocate buffer on aggregator node",
                    params->rank);
        }
    }

    status = MPI_Barrier(MPI_COMM_WORLD);

    comm_time += MPI_Wtime() - start_time;
    
    if (params->isWriter)
    {
        if (getenv("GCRM_LARGE_CHUNK")) {
            if (params->rank == 0) printf("using large chunk\n");
            if (params->nprocs % params->cb)
                ERROR("can't use large chunk with uneven CB decomposition",
                    params->rank);
            params->chunk[0] = params->nBuffered * params->ncells;
        } else if (getenv("GCRM_X_CHUNK")) {
            if (params->rank == 0) printf("using x-layer chunk\n");
            params->chunk[0] = params->nBuffered * params->ncells;
            params->chunk[1] = 1;
            params->chunk[2] = 1;
        }
        params->layout[1] = params->layout[0] + cbwidth * params->ncells - 1;

        if (params->chunking) {
            status = H5BlockDefine3DChunkDims (file,
                    params->chunk[0], params->chunk[1], params->chunk[2]);
            if (status != H5PART_SUCCESS) {
                ERROR("H5BlockDefine3DChunk failed", params->rank);
            }
        }

        status = H5BlockDefine3DFieldLayout (file,
                params->layout[0], params->layout[1],
                params->layout[2], params->layout[3],
                params->layout[4], params->layout[5]);
        if (status != H5PART_SUCCESS) {
            ERROR("H5BlockDefine3DFieldLayout failed", params->rank);
        }
    }

    for (i=0; i<params->t; i++)
    {
        // communication phase
        start_time = MPI_Wtime();

        status = MPI_Gather(
                data, (int)params->ndata, MPI_FLOAT,
                buffer, (int)params->ndata, MPI_FLOAT,
                0, comm);
        if (status != MPI_SUCCESS) ERROR("MPI_Gather failed", params->rank);

        comm_time += MPI_Wtime() - start_time;

        if (params->isWriter)
        {
            status = H5PartSetStep (file, i);
            if (status != H5PART_SUCCESS) {
                ERROR("H5PartSetStep failed", params->rank);
            }

            status = H5Block3dWriteScalarFieldFloat32 (file, "data", buffer);
            if (status != H5PART_SUCCESS) {
                ERROR("H5Block3dWriteScalarFieldFloat32 failed",
                        params->rank);
            }
        } 
    }

    if (params->isWriter) free(buffer);

    MPI_Group_free(&MPI_GROUP_WORLD);
    MPI_Group_free(&group);
    MPI_Comm_free(&comm);

    if (params->rank == 0) {
        printf ("(communication time = %g s)\n", comm_time);
    }
}

void write_h5part_cb3 (
        H5PartFile *file,
        float* data,
        Params* params)
{
    int i;
    double start_time, comm_time;
    h5part_int64_t status;
    int range[1][3];
    MPI_Group MPI_GROUP_WORLD;
    MPI_Group group;
    MPI_Comm comm;
    float *buffer;

    start_time = MPI_Wtime();
    comm_time = 0;

    int cbwidth = params->nBuffered;
    /*if ((params->myWriter + cbwidth) > params->nprocs)
                cbwidth = params->nprocs - params->rank;*/
    //printf("%d: writer %d cbwidth %d\n", params->rank, params->myWriter, cbwidth);

    // create new communicator for aggregation
    range[0][0] = params->myWriter; // start
    range[0][1] = params->myWriter + cbwidth - 1; // end
    range[0][2] = 1; // sride

    MPI_Comm_group(MPI_COMM_WORLD, &MPI_GROUP_WORLD);

    status = MPI_Group_range_incl(MPI_GROUP_WORLD, 1, range, &group);
    if (status != MPI_SUCCESS) ERROR("MPI_Group_range_incl failed", params->rank);

    status = MPI_Comm_create(MPI_COMM_WORLD, group, &comm);
    if (status != MPI_SUCCESS) ERROR("MPI_Comm_create failed", params->rank);

    if (params->isWriter) {
        buffer = (float *) malloc(params->nBuffered * params->ndata * sizeof(float));
        if (!buffer) {
            ERROR("Could not allocate buffer on aggregator node",
                    params->rank);
        }
    }

    status = MPI_Barrier(MPI_COMM_WORLD);

    comm_time += MPI_Wtime() - start_time;
       
    if (params->isWriter)
    {
        status = H5PartSetNumParticles (file,
                params->nBuffered * params->ndata);
        if (status != H5PART_SUCCESS) {
            ERROR("H5PartSetNumParticles failed", params->rank);
        }

        if (params->chunking) {
            status = H5PartSetChunkSize (file, params->ndata);
            if (status != H5PART_SUCCESS) {
                ERROR("H5PartSetChunkSize failed", params->rank);
            }
        }
    }

    for (i=0; i<params->t; i++)
    {
        // communication phase
        start_time = MPI_Wtime();

        status = MPI_Gather(
                data, (int)params->ndata, MPI_FLOAT,
                buffer, (int)params->ndata, MPI_FLOAT,
                0, comm);
        if (status != MPI_SUCCESS) ERROR("MPI_Gather failed", params->rank);

        comm_time += MPI_Wtime() - start_time;

        if (params->isWriter)
        {
            status = H5PartSetStep (file, i);
            if (status != H5PART_SUCCESS) {
                ERROR("H5PartSetStep failed", params->rank);
            }

            int j;
            for (j=0; j<params->nBuffered; j++)
            {
                h5part_int64_t offset = j*params->ndata;
                h5part_int64_t start = params->rank + offset;
                status = H5PartSetView (file, start, start + params->ndata - 1);
                if (status != H5PART_SUCCESS) {
                    ERROR("H5PartSetView failed", params->rank);
                }

                status = H5PartWriteDataFloat32 (file, "data", buffer + offset);
                if (status != H5PART_SUCCESS) {
                    ERROR("H5PartWriteDataFloat32 failed",
                            params->rank);
                }
            }
        } 
    }
    
    if (params->isWriter) free(buffer);

    MPI_Group_free(&MPI_GROUP_WORLD);
    MPI_Group_free(&group);
    MPI_Comm_free(&comm);

    if (params->rank == 0) {
        printf ("(communication time = %g s)\n", comm_time);
    }
}

double write_h5part (
        char* filename,
        float* data,
        Params* params)
{
    char flags;
    double start_time, end_time;
    H5PartFile* file;
    h5part_int64_t status;
    
    MPI_Barrier (MPI_COMM_WORLD);
    start_time = MPI_Wtime();
    
    flags = H5PART_WRITE;
    if (params->lustre) {
        flags |= H5PART_FS_LUSTRE;
    }
    if (params->mpiposix) {
        flags |= H5PART_VFD_MPIPOSIX;
    }

    if (params->isWriter)
    {
	/*
        file = H5PartOpenFileParallelAlign (
                filename,
                flags,
                params->commWriters,
                params->alignment);
	*/
        file = H5PartOpenFileParallel (
                filename,
                flags,
                params->commWriters);
        if (!file) {
            ERROR("Could not open H5Part file", params->rank);
        }

        if (params->throttle) {
            status = H5PartSetThrottle (file, params->throttle);
            if (status != H5PART_SUCCESS) {
                ERROR("H5PartSetThrottle failed", params->rank);
            }
        }
    }

    if (params->cb) {
        write_h5part_cb2 (file, data, params);
    } else {
        write_h5part_all (file, data, params);
    }

    if (params->isWriter) H5PartCloseFile (file);

    MPI_Barrier (MPI_COMM_WORLD);
    end_time = MPI_Wtime();

    return end_time - start_time;
}

void write_prs (Params* params)
{
    int i;
    int km;
    double time;
    double size;
    float* data;
    char* filename;

    filename = (char*) malloc (strlen(params->output) + 8);
    sprintf(filename, "%s/prs.h5", params->output);

    // the pressure variable has an additional "top" level
    km = params->km + 1;

    params->ncells  = params->im;
    params->ncells *= params->jm;
    params->ncells *= params->nm;

    params->ndata = params->ncells * km;

    data = (float*) malloc (params->ndata * sizeof(float));
    if (!data) {
        ERROR("Could not malloc data", params->rank);
    }

    for (i=0; i<params->ndata; i++) {
        data[i] = (float) random();
    }

    // each subdomain with its vertical levels is a chunk
    params->chunk[0] = params->im * params->jm;
    params->chunk[1] = km;
    params->chunk[2] = 1;

    params->layout[0] = params->rank * params->ncells;
    params->layout[1] = params->layout[0] + params->ncells - 1;
    params->layout[2] = 0;
    params->layout[3] = km - 1;
    params->layout[4] = 0;
    params->layout[5] = 0;

    if (params->rank == 0) remove (filename);


    time = write_h5part (filename, data, params);

    free (data);

    if (params->rank == 0) {
        size  = params->ngrid;
        size *= km;
        size *= sizeof(float);
        size /= ONE_MEGABYTE;
        size *= params->t;
        printf("Pressure variable (%d,%d,%d,%d)\n",
                params->im,
                params->jm,
                km,
                params->nm);
        printf("size = %g MB\n", size);
        printf("time = %g s\n", time);
        printf("bandwidth = %g MB/s\n", size / time);
    }
}

void write_edge (Params* params)
{
    int i;
    double time;
    double size;
    float* data;
    char* filename;

    filename = (char*) malloc (strlen(params->output) + 8);
    sprintf(filename, "%s/edge.h5", params->output);

    params->ncells  = params->im;
    params->ncells *= params->jm;
    params->ncells *= params->nm;

    params->ndata = 6 * params->ncells * params->km;

    data = (float*) malloc (params->ndata * sizeof(float));
    if (!data) {
        ERROR("Could not malloc data", params->rank);
    }

    for (i=0; i<params->ndata; i++) {
        data[i] = (float) random();
    }

    // each subdomain with its vertical levels is a chunk
    params->chunk[0] = params->im * params->jm;
    params->chunk[1] = params->km;
    params->chunk[2] = 6;

    params->layout[0] = params->rank * params->ncells;
    params->layout[1] = params->layout[0] + params->ncells - 1;
    params->layout[2] = 0;
    params->layout[3] = params->km - 1;
    params->layout[4] = 0;
    params->layout[5] = 5;

    if (params->rank == 0) remove (filename);

    time = write_h5part (filename, data, params);

    free (data);

    if (params->rank == 0) {
        size  = params->ngrid;
        size *= params->km;
        size *= sizeof(float);
        size *= 6.0;
        size /= ONE_MEGABYTE;
        size *= params->t;
        printf("Edge-centered variable (%d,%d,%d,%d,6)\n",
                params->im,
                params->jm,
                params->km,
                params->nm);
        printf("size = %g MB\n", size);
        printf("time = %g s\n", time);
        printf("bandwidth = %g MB/s\n", size / time);
    }
}


/******************************************************************************/
/* Main procedure                                                             */
/******************************************************************************/

int main (int argc, char** argv)
{
    Params params;

    // initialize MPI
    MPI_Init (&argc, &argv);
    MPI_Comm_rank (MPI_COMM_WORLD, &params.rank);
    MPI_Comm_size (MPI_COMM_WORLD, &params.nprocs);

    if (params.rank == 0) {
        printf ("H5_VERS_INFO: %s\n", H5_VERS_INFO);
        printf ("H5PART_VER_STRING: %s\n", H5PART_VER_STRING);
    }

    parse_args (argc, argv, &params);

    H5PartSetVerbosityLevel (params.verbosity);

    if (params.cb) {
        init_cb2(&params);
    } else {
        params.isWriter = 1;
        params.commWriters = MPI_COMM_WORLD;
    }

    int i;
    for (i=0; i<params.i; i++)
    {
        if (params.rank == 0) {
            printf ("Write test %d...\n", i);
        }

        if (params.prs)  write_prs  (&params);
        if (params.edge) write_edge (&params);
        
        if (!(params.prs) && !(params.edge)) {
            if (params.rank == 0) {
                printf ("No variables specified... exiting without write!\n");
            }
            break;
        }
    }
    
    MPI_Finalize();
    return (EXIT_SUCCESS);
}

