dnl /**
dnl *
dnl *** Copyright Notice ***
dnl
dnl PIOK - Parallel I/O Kernels, Copyright (c) 2015-2016, 
dnl The Regents of the University of California, through Lawrence
dnl Berkeley National Laboratory (subject to receipt of any required approvals
dnl from the U.S. Dept. of Energy).  All rights reserved.
dnl  
dnl If you have questions about your rights to use or distribute this software,
dnl please contact Berkeley Lab's Innovation & Partnerships Office at
dnl IPO@lbl.gov.
dnl  
dnl NOTICE.  This Software was developed under funding from the U.S. Department
dnl of Energy and the U.S. Government consequently retains certain rights. As
dnl such, the U.S. Government has been granted for itself and others acting on
dnl its behalf a paid-up, nonexclusive, irrevocable, worldwide license in the
dnl Software to reproduce, distribute copies to the public, prepare derivative
dnl works, and perform publicly and display publicly, and to permit other to do
dnl so.
dnl  
dnl */
dnl /**
dnl *
dnl * Email questions to SByna@lbl.gov
dnl * Scientific Data Management Research Group
dnl * Lawrence Berkeley National Laboratory
dnl *
dnl * Fri Mar 25 09:59:39 PDT 2016
dnl *
dnl */

AC_DEFUN(my_HDF5_PATH, [
HDF5_PATH=
AC_ARG_WITH(hdf5,
[  --with-hdf5=<PATH>          Specify the HDF5 path],
if test -d "$withval" -a -d $withval/include/ -a -d $withval/lib ; then
    if test -f $withval/include/hdf5.h -a -f $withval/bin/h5pcc ; then
        HDF5_PATH=$withval
    else
        AC_MSG_ERROR([invalid HDF5 path specified])
    fi
fi
)
if test "$HDF5_PATH" = "" ; then
    AC_MSG_CHECKING([for right HDF5 path])
    tmp=`/usr/bin/which h5dump 2>&1`
    if test -f "$tmp" ; then
        tmp=`echo $tmp | sed 's%/bin/h5dump%%'`
        if test -d $tmp/include -a -d $tmp/lib ; then
            if test -f $tmp/include/hdf5.h ; then
                HDF5_PATH=$tmp
                AC_MSG_RESULT([HDF5 found])
            fi
        fi
    else
        AC_MSG_RESULT([HDF5 not found])
    fi
fi
if test "$HDF5_PATH" = "" ; then
    AC_MSG_ERROR([invalid HDF5 PATH])
fi
])

AC_DEFUN(my_H5PART_PATH, [
H5PART_PATH=
AC_ARG_WITH(h5part,
[  --with-h5part=<PATH>          Specify the H5Part path],
if test -d "$withval" -a -d $withval/include/ -a -d $withval/lib ; then
    if test -f $withval/include/H5Part.h ; then
        H5PART_PATH=$withval
    else        
        AC_MSG_ERROR([invalid H5Part path specified])
    fi
fi
)
if test "$H5PART_PATH" = "" ; then
    AC_MSG_ERROR([invalid H5PART PATH])
fi
])

AC_DEFUN(my_MPI_PATH, [
MPI_PATH=
AC_ARG_WITH(mpi,
[  --with-mpi=<PATH>          Specify the MPI path],
if test -d "$withval" -a -d $withval/include/ -a -d $withval/lib ; then
    if test -f $withval/include/mpi.h ; then
        MPI_PATH=$withval
    else
        AC_MSG_ERROR([invalid MPI path specified])
    fi
fi
)
if test "$MPI_PATH" = "" ; then
    AC_MSG_CHECKING([for right MPI path])
    tmp=`/usr/bin/which mpicc 2>&1`
    if test -f "$tmp" ; then
        tmp=`echo $tmp | sed 's%/bin/mpicc%%'`
        if test -d $tmp/include -a -d $tmp/lib ; then
            if test -f $tmp/include/mpi.h ; then
                MPI_PATH=$tmp
                AC_MSG_RESULT([MPI found])
            fi
        fi
    else
        AC_MSG_RESULT([MPI not found])
    fi
fi
if test "$MPI_PATH" = "" ; then
    AC_MSG_WARN([invalid MPI PATH. will use system default if any])
fi
])

AC_DEFUN(my_GET_INSTALL, [
INSTALL=
tmp=`which install 2>&1`
AC_MSG_CHECKING([for install])
if test -f "$tmp" ; then
	INSTALL=$tmp
	AC_MSG_RESULT([ok])
else
	AC_MSG_WARN([install does not exist])
	INSTALL="cp -f"
	AC_MSG_RESULT([ok])
fi
])

AC_DEFUN(my_CC_GET_COMPILER, [
CC_COMPILER=
tmp=`which cc 2>&1`

AC_MSG_CHECKING([for Cray cc compiler])
if test -f $tmp ; then
    CC_COMPILER=$tmp
    AC_MSG_RESULT([found])
elif test -n $CC ; then
    CC_COMPILER=$CC
    AC_MSG_RESULT([found, but not necessarily Cray cc])
else
    AC_MSG_RESULT([not found])
    AC_MSG_ERROR([valid Cray cc compiler required])
fi
])

AC_DEFUN(my_CCC_GET_COMPILER, [
CCC_COMPILER=
tmp=`which CC 2>&1`

AC_MSG_CHECKING([for Cray CC compiler])
if test -f $tmp ; then
    CCC_COMPILER=$tmp
    AC_MSG_RESULT([found])
elif test -n $CC ; then
    CCC_COMPILER=$CC
    AC_MSG_RESULT([found, but not necessarily Cray CC])
else
    AC_MSG_RESULT([not found])
    AC_MSG_ERROR([valid Cray CC compiler required])
fi
])

AC_DEFUN(my_H5PCC_GET_COMPILER, [
H5PCC_COMPILER=
tmp=`which h5pcc 2>&1`

AC_MSG_CHECKING([for H5PCC compiler])
if test -f $tmp ; then
    H5PCC_COMPILER=$tmp
    AC_MSG_RESULT([found])
elif test -n $H5PCC ; then
    H5PCC_COMPILER=$H5PCC
    AC_MSG_RESULT([found, but not necessarily H5PCC])
else
    AC_MSG_RESULT([not found])
    AC_MSG_ERROR([valid H5PCC compiler required])
fi
])


AC_DEFUN(my_GCC_GET_COMPILER, [
GCC_COMPILER=
tmp=`which gcc 2>&1`

AC_MSG_CHECKING([for GNU C compiler])
if test -f $tmp ; then
	GCC_COMPILER=$tmp
	AC_MSG_RESULT([found])
elif test -n $CC ; then
    GCC_COMPILER=$CC
    AC_MSG_RESULT([found, but not necessarily gcc])
else
	AC_MSG_RESULT([not found])
	AC_MSG_ERROR([valid gcc compiler required])
fi
])

AC_DEFUN(my_GXX_GET_COMPILER, [
GXX_COMPILER=
tmp=`which g++ 2>&1`

AC_MSG_CHECKING([for GNU C++ compiler])
if test -f $tmp ; then
    GXX_COMPILER=$tmp
    AC_MSG_RESULT([found])
elif test -n $CXX ; then
    GXX_COMPILER=$CXX
    AC_MSG_RESULT([found, but not necessarily g++])
else
    AC_MSG_RESULT([valid g++ compiler not found])
fi
])

