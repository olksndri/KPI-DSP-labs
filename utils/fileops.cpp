#include "fileops.h"

void log_to_file(char *logname, void *data, int len, DATA_TYPE type) 
{   
    char fname[512]; 
    snprintf(fname, 512, "logs/%s.txt", logname); 

    FILE *f_log = fopen(fname, "w"); 

    if(f_log == NULL)
    { 
        perror("Error opening file! Aborting\n");
        abort();  
    }   

    if(type == COMPLEX) 
    {
        complex *p_data = (complex*)data; 
        for(int i = 0; i < len; i++)
            fprintf(f_log, "%d  real: %f \t imag: %f\n", i, p_data[i].real, p_data[i].imag); 
    }
    else 
    {  
        perror("Not implemented! Aborting\n");
        abort();  
    }

    fclose(f_log); 

}