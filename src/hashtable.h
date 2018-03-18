#ifndef _HASHTABLE_
#define _HASHTABLE_

#include "functions.h"

// Declare functions of hashtable.c file
hashtable **hashtable_init(int, int);

void hashtable_insert(hashtable *, int, node *);

void hashtable_destroy(hashtable **, int);

int hash_function_classic(double *, int, int, int *);

int **generate_array_r_classic(int, int);

int hash_function_LSH(double *, int, int, int *, int, double **);

double **generate_array_LSH(int, int);

void hashtable_print(hashtable *);

#endif