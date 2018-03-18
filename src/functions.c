#include "functions.h"

// Function that creats the grid curve of a curve
double** find_grid_curve(curve *real_curve, double d, double t)
{
	int i, j;
	double a;
	double **grid_curve;
		
	// Creation of an array that is practically the grid curve , with the same size as the real curve
	grid_curve = (double**) malloc((real_curve->points) * sizeof(double*));
	if (grid_curve == NULL)
	{
		printf("Malloc: memory allocation error!\n");
		exit(3);
	}
	
	for (i = 0 ; i < (real_curve->points) ; i++)
	{
		grid_curve[i] = (double*) malloc((real_curve->dimension) * sizeof(double));
		
		if (grid_curve[i] == NULL)
		{
			printf("Malloc: memory allocation error!\n");
			exit(3);
		}
	}
			
	for(i = 0; i < (real_curve->points); i++)
	{
		for(j = 0; j < (real_curve->dimension); j++)
		{
			a = round(real_curve->coordinates[i][j] / d);
						
			grid_curve[i][j] = (double)((int) a * d) + t;
		}	
	}
	return grid_curve;
}


/* Function that removes consecutive duplicate coordinates , transforms multi-dimensional grid curve to one-dimensional 
grid curve and creates a new node that consists of: 
					1) info about initial curve
					2) size of new grid curve ( that has no more duplicates )
					3) the grid curve that has no more duplicates
					4) the pointer to the next node */
node *create_hashtable_node(double **grid_curve, curve *real_curve, int K) 
{
    int i, j, size = 1, distinct;
    double **grid_no_duplicate;
	double *grid_curve_1D;

	// Malloc for the first line of our new array
    grid_no_duplicate = (double**) malloc(size * sizeof(double*));
	
	if (grid_no_duplicate == NULL)
	{
		printf("Malloc: memory allocation error!\n");
		exit(3);
	}
		
	// Malloc for the number of dimensions of our vector
    grid_no_duplicate[0] = (double*) malloc((real_curve->dimension) * sizeof(double));
	
	if (grid_no_duplicate[0] == NULL)
	{
		printf("Malloc: memory allocation error!\n");
		exit(3);
	}
    
    grid_no_duplicate[0] = grid_curve[0]; 

	// Comparison between the coordinates of points of our old grid_curve
    for (i = 1 ; i < (K * (real_curve->points)) ; i++)
    {
        distinct = 0;
		
		// If there is at least one different coordinate, the vectors will be different
        for (j = 0 ; j < (real_curve->dimension) ; j++)
        {
            if (grid_curve[i][j] !=  grid_curve[i-1][j])
            {
				distinct = 1;
				break;
			}
        }

		// If the coordinates weren't identically the same , we realloc our new array to append the new point
        if (distinct)
        {
            grid_no_duplicate = (double**) realloc(grid_no_duplicate, (++size) * sizeof(double*));
			
			if (grid_no_duplicate == NULL)
			{
				printf("Realloc: memory allocation error!\n");
				exit(3);
			}
			
            grid_no_duplicate[size-1] = (double*) malloc((real_curve->dimension)  * sizeof(double));
			
			if (grid_no_duplicate[size-1] == NULL)
			{
				printf("Malloc: memory allocation error!\n");
				exit(3);
			}
            
            grid_no_duplicate[size-1] = grid_curve[i];   
        }
    }
	
	// Transform grid curve to one dimension vector
	grid_curve_1D = grid_curve_to_1D(grid_no_duplicate, size, (real_curve->dimension));
	
	
	// Free space for grid_no_duplicate
	for (i = 0 ; i < size ; i++)
		free(grid_no_duplicate[i]);
	free(grid_no_duplicate);
	
	
	// New node for hashtable's linked list
	node *new_node;
	
	// Allocate memory for new_node
	new_node = (node*) malloc(sizeof(node));
	
	if (new_node == NULL)
	{
		printf("Malloc: memory allocation error!\n");
		exit(3);
	}
	
	// Fill struct with data
	new_node->real_curve = real_curve;
	new_node->grid_points = (size * (real_curve->dimension));
	new_node->grid_curve_1D = grid_curve_1D;
	new_node->next = NULL;
	
	return new_node;
}


/* Function that transforms 2D array of grid curve to 1D array 
			- points: number of grid curve points after duplicates removal */
double *grid_curve_to_1D(double **grid_curve, int points, int dimension)
{
	int i, j;
	double *grid_curve_1D;
	
	grid_curve_1D = (double*) malloc((points * dimension) * sizeof(double));
	
	if (grid_curve_1D == NULL)
	{
		printf("Malloc: memory allocation error!\n");
		exit(3);
	}
	
	for (i = 0 ; i < points ; i++)
		for (j = 0 ; j < dimension ; j++)			
			grid_curve_1D[(i * dimension) + j] = grid_curve[i][j];
	
	return grid_curve_1D;
}