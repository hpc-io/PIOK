PIOK: Parallel I/O Kernels 
==========================
Suren Byna and Mark Howison
Contact: SByna@lbl.gov

Parallel I/O is an essential part of many scientific simulation and
analysis codes on high performance computing (HPC) systems, as the codes
need to store or retrieve data to or from parallel file systems,
respectively. Identifying I/O performance bottlenecks, stressing parallel
file systems, and tuning performance of real scientific applications is
critical in using HPC systems efficiently. Instead of running large-scale
simulations that contain full-blown computation and communication phases
along with I/O, extracting the parallel I/O portions of codes reduces the
complexity in tuning I/O. Toward that goal, we extract I/O kernels from
real scientific simulations and call the set of kernels “Parallel I/O
Kernels (PIOK)”. As of September 2015, we have developed three I/O kernels
from three scientific simulations: VPIC-IO, VORPAL-IO, and GCRM-IO. 

VPIC-IO (1D array): VPIC is a highly optimized and scalable particle
physics simulation developed by Los Alamos National Lab [1]. VPIC-IO uses
the H5Part API to create a file, write eight variables, and close the file.
H5Part provides a simple veneer API for issuing HDF5 calls corresponding to
a time-varying, multivariate particle data model [2]. We extracted all the
H5Part function calls of the VPIC code to form the VPIC-IO kernel. The
particle data written in the kernel is random data of oat data type. The
I/O pattern of VPIC-IO is a 1D particle array of a given number of
particles and each particle has eight variables. The kernel, by default,
writes 8M particles per MPI process for all experiments reported in this
paper.

VORPAL-IO (3D block structured grid): This I/O kernel is extracted from
VORPAL, a computational plasma framework application simulating the
dynamics of electromagnetic systems, plasmas, and rarefied as well as dense
gases, developed by TechX [3]. This benchmark uses H5Block to write
non-uniform chunks of 3D data per processor. The kernel takes 3D block
dimensions (x, y, and z) and the number of components as input.

GCRM-IO (semi structured mesh): This I/O kernel simulates I/O for GCRM, a
global atmospheric circulation model, simulating the circulations
associated with large convective clouds. This I/O kernel also uses H5Part
to perform I/O operations. The kernel performs all the GCRM I/O operations
with random data. I/O pattern of GCRM-IO corresponds to a semi-structured
geodesic mesh, where the grid resolution and sub-domain resolution are
specified as input. 

[1]	E. W. Bethel, J. M. Shalf, C. Siegerist, K. Stockinger, A.
Adelmann, A. Gsell, B. Oswald, and T. Schietinger. Progress on H5Part:  A
Portable High Performance Parallel Data Interface for Electromagnetics
Simulations. In Proceedings of the 2007 IEEE Particle Accelerator
Conference (PAC 07). 25-29 Jun 2007, Albuquerque, New Mexico. 22nd IEEE
Particle Accelerator Conference, p.3396 , 2007.
[2]	K. J. Bowers, B. J. Albright, L. Yin, B. Bergen, and T. J. T. Kwan.
Ultrahigh performance three-dimensional electromagnetic relativistic
kinetic plasma simulation. Physics of Plasmas, 15(5):7, 2008.
[3]	C. Nieter and J. R. Cary. VORPAL: a versatile plasma simulation
code. Journal of Computational Physics, 196:448{472, 2004.
[4]	David A. Randall et al., Global cloud resolving model (GCRM),
http://kiwi.atmos.colostate.edu/gcrm/ 

-------------------------------------------------------------------------

REQUIREMENTS
=============

This code has been tested with gcc on NERSC's Edison and Hopper
supercomputers. We have used gcc compiler to build them. Compiling the
kernels need HDF5 and H5Part libraries. We have compiled HDF5 with
--enable-parallel flag.

Example for compiling HDF5 1.8.14 and H5Part on Edison:
> cd <HDF5_source_dir>
> CC=cc ./configure --prefix=<HDF5_Install_dir> --disable-shared --disable-production --enable-parallel
> make
> make install

> cd <H5Part_source_dir>
> CC=cc ./configure --prefix=$PWD --with-hdf5=<HDF5_Install_dir> --enable-parallel --no-create --no-recursion
> make
> make install

-------------------------------------------------------------------------

CODE STRUCTURE
===============

The kernels are arranged in different directories. Currently, we have
gcrmio, vorpalio, and vpicio directories.

In <piok_dir>/gcrmio, we have the code for gcrmio kernel. 
<piok_dir>/gcrmio/README file has further details on compiling, and sample
commands to run on NERSC's Edison system.

In <piok_dir>/vorpalio, we have the code for vorpalio kernel. 
<piok_dir>/vorpalio/README file has further details on compiling, and
sample commands to run on NERSC's Edison system.

In <piok_dir>/vpicio, we have the code for multiple versions of vpicio kernel. 
	<piok_dir>/vpicio/vpicio_hdf5: Uses HDF5 calls to perform file
				       writes.
	<piok_dir>/vpicio/vpicio_uni:  Uses H5Part calls to perform file
				       writes, which internally calls HDF5.
				       _uni refers to each MPI process
					writing	same number of particles.
	<piok_dir>/vpicio/vpicio_non_uni:  Uses H5Part calls to perform file
				       writes, which internally calls HDF5.
				       _non_uni refers to each MPI process
					writing	different number of
					particles. The variance in the
					number of particles is controlled
					by VARIABILITY constant.
	<piok_dir>/vpicio/vpicio_uni_md:  Uses HDF5's multi-dataset write
					functions. Multi-dataset functions
					allow writing multiple datasets using
					a single write call. 

	Compiling and running vpicio kernels are simple. Modify the HDF5
	and H5Part paths in the corresponding Makefile and run "make". 
	<piok_dir>/vpicio/README file has further details on compiling, and sample
	commands to run on NERSC's Edison system.

