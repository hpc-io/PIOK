#ifndef READH5_H
#define READH5_H

#include "params.h"

long int read_vapor_data_nc_1day (float * vapor_data, input_params * ips, \
                                  int rows_from, int rows_to,           \
                                  int cols_from, int cols_to);


#endif

