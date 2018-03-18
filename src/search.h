#ifndef _SEARCH_
#define _SEARCH_

#include "functions.h"

// Declare functions of search.c file
void search(char *, char *, database *, int, int, int, int);

void search_stats(char *, database *, double **, int, int, int, int);

int compare_grid_curves(double *, double *, int);

tuple_true_neighbor *find_true_nearest_neighbor(hashtable *, node *, int, int);

int *remove_duplicates(int *, int);

int number_of_curves_query(char *);

#endif