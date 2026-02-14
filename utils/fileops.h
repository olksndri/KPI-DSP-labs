#pragma once

#include <stdlib.h>
#include <stdio.h> 

#include "../math/complex_math.h" 

#include "utils.h" 
// #include <fcntl.h> 

void log_to_file(char *logname, void *data, int len, DATA_TYPE type); 
