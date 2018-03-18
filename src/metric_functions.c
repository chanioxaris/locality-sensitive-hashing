#include "functions.h"

// Function that finds out the max value between two numbers
double max_2(double a, double b)
{
	if (a >= b)
    	return a;
  
  	return b;
}


// Function that finds out the min value out of three numbers
double min_3(double a, double b, double c)
{
	if (a <= b && a <= c)
      	return a;
  	else if (b <= a && b <= c)
      	return b;
    else
      	return c;
}


// Function that calculates the Euclidean distance between two points
double euclidean(double *coord1, double *coord2, int dimension)
{
	int i;
  	double dist = 0.0;
  
  	for (i = 0 ; i < dimension ; i++)
      	dist += pow((coord1[i] - coord2[i]), 2);
  
	return sqrt(dist);
}


// Function that calculates the Frechet distance between two curves (using a 2D array)
double frechet_distance(double **curve1, double **curve2, int m1, int m2, int dimension)
{
  	int i, j;
	double distance;
	double **array;
  	
  	array = (double**) malloc(m1 * sizeof(double*));
	
	if (array == NULL)
	{
		printf("Malloc: memory allocation error!\n");
		exit(3);
	}
 
    for (i = 0 ; i < m1 ; i++)
	{
    	array[i] = (double*) malloc(m2 * sizeof(double));
		
		if (array[i] == NULL)
		{
			printf("Malloc: memory allocation error!\n");
			exit(3);
		}
	}

    array[0][0] = euclidean(curve1[0], curve2[0], dimension);
	
    for (i = 1 ; i < m1 ; i++)
    	array[i][0] = max_2(array[i-1][0], euclidean(curve1[i], curve2[0], dimension));
	
    for (j = 1 ; j < m2 ; j++)
    	array[0][j] = max_2(array[0][j-1], euclidean(curve1[0], curve2[j], dimension));
	
    for (i = 1 ; i < m1 ; i++)
      	for (j = 1 ; j < m2 ; j++)
        	array[i][j] = max_2(min_3(array[i-1][j], array[i][j-1], array[i-1][j-1]), euclidean(curve1[i], curve2[j], dimension));
	
	distance = array[m1-1][m2-1];
	
    for (i = 0 ; i < m1 ; i++)
    	free(array[i]);      
    free(array);
	
	return distance;
}


// Function that calculates the DTW distance between two curves (using a 2D array)
double DTW_distance(double **curve1, double **curve2, int m1, int m2, int dimension)
{
	int i, j;
	double distance;
	double **array;
  		
  	array = (double**) malloc((m1+1) * sizeof(double*));
	
	if (array == NULL)
	{
		printf("Malloc: memory allocation error!\n");
		exit(3);
	}
      
    for (i = 0 ; i <= m1 ; i++)
	{
    	array[i] = (double*) malloc((m2+1) * sizeof(double));
		
		if (array[i] == NULL)
		{
			printf("Malloc: memory allocation error!\n");
			exit(3);
		}		
	}

  	array[0][0] = 0.0;
  	
  	for (i = 1 ; i <= m1 ; i++)
    	array[i][0] = INFINITY;
      	  
    for (j = 1 ; j <= m2 ; j++)
    	array[0][j] = INFINITY;
 
  	for (i = 1 ; i <= m1 ; i++)
      	for (j = 1 ; j <= m2 ; j++)
        	array[i][j] = euclidean(curve1[i-1], curve2[j-1], dimension) + min_3(array[i-1][j], array[i][j-1], array[i-1][j-1]);
  

	distance = array[m1][m2];
  
  	for (i = 0 ; i <= m1 ; i++)
    	free(array[i]);      
    free(array);
      
	return distance;
}	