#include "memops.h"

void* malloc_nc(size_t size) 
{
    void *ptr = malloc(size);  

	if(ptr == NULL)
	{ 
		printf("Error in memory allocation! Aborting program... \n"); 
		abort();
	} 

	return ptr; 
}

void free_nc(void *ptr)
{
	if(ptr != NULL)
		free(ptr); 
}