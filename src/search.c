#include "functions.h"

double frechet_distance(double **, double **, int, int, int);

double DTW_distance(double **, double **, int, int, int);


// Function that searches for query curve 
void search(char *query_file, char *output_file, database *data_search, int L, int K, int hash_function, int metric_function)
{
  	int i, j, z, index, min_ID, found_grid, size_neigh_ID;
	int *neigh_ID, *indexes, *neigh_ID_no_duplicates;
	double R, distance, min_distance;
	char line[LINE_LENGTH];
  	char *data;
	double **tmp_grid_curve_conc, **tmp_grid_curve;
	node *new_node, *tmp;
	curve *tmp_real;
	tuple_true_neighbor *true_neigh;

	// Open input file
	FILE *query = fopen(query_file, "r");
  
  	if (query == NULL)
  	{
    	printf("Fopen: error opening query file!\n");
      	exit(1);
    }
		
	// Get diameter from query file
  	fgets(line, sizeof(line), query);
  	strtok(line,"\t");
  	
	data = strtok(NULL,"\n "); 
  	R = atof(data);
	
	// Allocate memory for array to store buckets indexes
	indexes = (int*) malloc(L * sizeof(int));
	
	if (indexes == NULL)
	{
		printf("Malloc: memory allocation error!\n");
		exit(3);
	}

  	while (fgets(line, sizeof(line), query) != NULL)
    {			
		curve *new_curve = (curve*) malloc(sizeof(curve));
		
		if (new_curve == NULL)
		{
			printf("Malloc: memory allocation error!\n");
			exit(3);
		}
		
		// Assign curve dimension
		new_curve->dimension = data_search->dimension;
		
    	data = strtok(line,"\t");
		
		// Assign curve's ID
		new_curve->ID_curve = atoi(data);

		data = strtok(NULL,"\t");
		
		// Assign the points of curve
		new_curve->points = atoi(data);
		
		// Allocate memory for array to store the points
		new_curve->coordinates = (double**) malloc((new_curve->points) * sizeof(double*));
		
		if (new_curve->coordinates == NULL)
		{
			printf("Malloc: memory allocation error!\n");
			exit(3);
		}

		// Parse every single point of curve
      	for (i = 0 ; i < new_curve->points ; i++)
        {	
			// Allocate memory for array to store each coordinate for every point (related to dimension {2,3,4})	
			new_curve->coordinates[i] = (double*) malloc(data_search->dimension * sizeof(double));
			
			if (new_curve->coordinates[i] == NULL)
			{
				printf("Malloc: memory allocation error!\n");
				exit(3);
			}
			
			// Cases for different dimensions
			switch(data_search->dimension)
			{
				case 2:
					data = strtok(NULL,",");
				
					// Assign 1st Coordinate
					data = data + 1;				
					new_curve->coordinates[i][0] = atof(data);
					
					// Assign 2st Coordinate
					data = strtok(NULL,")");				
					new_curve->coordinates[i][1] = atof(data);

					strtok(NULL," ");		
					break;
				case 3:
					data = strtok(NULL,",");			
				
					// Assign 1st Coordinate
					data = data + 1;
					new_curve->coordinates[i][0] = atof(data);

					// Assign 2st Coordinate
					data = strtok(NULL,",");
					new_curve->coordinates[i][1] = atof(data);
					
					// Assign 3st Coordinate
					data = strtok(NULL,")");
					new_curve->coordinates[i][2] = atof(data);				
					
					strtok(NULL," ");			
					break;
				case 4:
					data = strtok(NULL,",");
				
					// Assign 1st Coordinate
					data = data + 1;
					new_curve->coordinates[i][0] = atof(data);
				  
					// Assign 2st Coordinate
					data = strtok(NULL,",");
					new_curve->coordinates[i][1] = atof(data);
				  
					// Assign 3st Coordinate
					data = strtok(NULL,",");
					new_curve->coordinates[i][2] = atof(data);
					
					// Assign 4st Coordinate
					data = strtok(NULL,")");
					new_curve->coordinates[i][3] = atof(data);
					
					strtok(NULL," ");			
					break;
				default:
					printf("Wrond input dimension!\n");
					exit(1);	
			}	
        }

		// Temporary array to concatenate all grid curves
		tmp_grid_curve_conc = (double**) malloc((K * (new_curve->points)) * sizeof(double*));
		
		if (tmp_grid_curve_conc == NULL)
		{
			printf("Malloc: memory allocation error!\n");
			exit(3);
		}
		
		for (i = 0 ; i < (K * (new_curve->points)) ; i++)
		{
			tmp_grid_curve_conc[i] = (double*) malloc((new_curve->dimension) * sizeof(double));
			
			if (tmp_grid_curve_conc[i] == NULL)
			{
				printf("Malloc: memory allocation error!\n");
				exit(3);
			}
		}
		
		size_neigh_ID = 0;
		found_grid = 0;
		min_distance = INFINITY;
		min_ID = -1;
		
		neigh_ID = (int*) malloc(size_neigh_ID * sizeof(int));

		// Fill L hashtables
		for (i = 0 ; i < L ; i++)
		{				
			// Find k grid curves
			for (j = 0 ; j < K ; j++)
			{			  
				// Find grid curve for each grid
				tmp_grid_curve = find_grid_curve(new_curve, data_search->d, data_search->array_grid_t[i][j]);
				
				for (z = (j * (new_curve->points)) ; z < ((j+1) * (new_curve->points)) ; z++)
					tmp_grid_curve_conc[z] = tmp_grid_curve[z - (j * (new_curve->points))];
			}
			
			// Create final grid vector	and store to struct _node_		
			new_node = create_hashtable_node(tmp_grid_curve_conc, new_curve, K);
			
			
			// Check which hash funtion to use (CLASSIC or LSH)
			if (hash_function == CLASSIC)
				index = hash_function_classic(new_node->grid_curve_1D, new_node->grid_points, 
											data_search->ht[i]->size, data_search->array_r_classic[i]);				
			else
				index = hash_function_LSH(new_node->grid_curve_1D, new_node->grid_points, data_search->ht[i]->size, 
											data_search->array_r_classic[i], i, data_search->array_LSH);

			// Fill the bucket's index for L's hashtable
			indexes[i] = index;
			
			tmp = data_search->ht[i]->table[index];
			
			// Search the whole list of this table[index]
			while (tmp != NULL)
			{
				// Only for curves with the same amount of grid points
				if (new_node->grid_points == tmp->grid_points)
				{	
					// Only for curves with same 1D grid curves
					if (compare_grid_curves(new_node->grid_curve_1D, tmp->grid_curve_1D, tmp->grid_points))
					{
						
						found_grid = 1;
						
						// Calculate distance
						if (metric_function == DFT)
							distance = frechet_distance(new_node->real_curve->coordinates, tmp->real_curve->coordinates, 
										new_node->real_curve->points, tmp->real_curve->points, tmp->real_curve->dimension);
						else
							distance = DTW_distance(new_node->real_curve->coordinates, tmp->real_curve->coordinates, 
										new_node->real_curve->points, tmp->real_curve->points, tmp->real_curve->dimension);	

						// Nearest neighbor
						if (distance < min_distance)
						{
							min_distance = distance;
							min_ID = tmp->real_curve->ID_curve;
						}

						// R diameter neighbors						
						if (distance < R)
						{
							// Every single time that there is a "near-neighbor" ,memory is reallocated..
							neigh_ID = (int*) realloc(neigh_ID, (++size_neigh_ID) * sizeof(int));
							
							if (neigh_ID == NULL)
							{
								printf("Realloc: memory allocation error1!\n");
								exit(3);
							}
							
							//.. to store it
							neigh_ID[size_neigh_ID-1] = tmp->real_curve->ID_curve;
						}
					}
				}
				// Step to the next node in bucket
				tmp = (tmp->next);
			}	
		}	
				
		// If we dont find any same grid curve in buckets, check all the curves in these buckets
		if (!found_grid)
		{
			// For L hashtables
			for (i = 0 ; i < L ; i++)
			{
				// Assign for the i' hashtable the i' index of query curve
				tmp = data_search->ht[i]->table[indexes[i]];
				
				// For the whole bucket with this index
				while (tmp != NULL)
				{
					// Calculate distance
					if (metric_function == DFT)
						distance = frechet_distance(new_node->real_curve->coordinates, tmp->real_curve->coordinates, 
									new_node->real_curve->points, tmp->real_curve->points, tmp->real_curve->dimension);
					else
						distance = DTW_distance(new_node->real_curve->coordinates, tmp->real_curve->coordinates, 
									new_node->real_curve->points, tmp->real_curve->points, tmp->real_curve->dimension);
					
					
					// Nearest neighbor
					if (distance < min_distance)
					{
						min_distance = distance;
						min_ID = tmp->real_curve->ID_curve;
					}
												
					// R diameter neighbors						
					if (distance < R)
					{
						neigh_ID = (int*) realloc(neigh_ID, (++size_neigh_ID) * sizeof(int));
						
						if (neigh_ID == NULL)
						{
							printf("Realloc: memory allocation error2!\n");
							exit(3);
						}
						
						neigh_ID[size_neigh_ID-1] = tmp->real_curve->ID_curve;
					}
					
					// Step
					tmp = (tmp->next);
				}
			}
		}
			
		// Sort the array that contains R-neighbor's ID 
		quickSort(neigh_ID, 0, size_neigh_ID-1);
		
		// Remove duplicates from neigh_ID
		neigh_ID_no_duplicates = remove_duplicates(neigh_ID, size_neigh_ID);
		
		// Find the true neighbor for the curve
		true_neigh = find_true_nearest_neighbor(data_search->ht[0], new_node, L, metric_function);
		
		//OUTPUT FUNCTION
		output_nonstats(output_file, new_node->real_curve->ID_curve, metric_function, hash_function, found_grid,
                     min_ID, true_neigh->true_min_ID, min_distance, true_neigh->true_min_distance, neigh_ID_no_duplicates);	
				
				
		// Free new_node			
		tmp_real = new_node->real_curve;

		for (z = 0 ; z < tmp_real->points ; z++)
			free(tmp_real->coordinates[z]);
		free(tmp_real->coordinates);

		free(tmp_real);	
		
		free(new_node->grid_curve_1D);
		free(new_node);
		
		// Free neigh_ID_no_duplicates	
		free(neigh_ID_no_duplicates);
				
		// Free true_neigh
		free(true_neigh);
	}
		
	// Free indexes
	free(indexes);
		
	fclose(query);
	return;
}


// Function that searches for query curve with parameter -stats
void search_stats(char *query_file, database *data_search, double **stats, int L, int K, int hash_function, int metric_function)
{
  	int i, j, z, index, counter = 0, min_ID, found_grid, size_neigh_ID;
	int *neigh_ID, *indexes, *neigh_ID_no_duplicates;
	double R, distance, min_distance;
	char line[LINE_LENGTH];
  	char *data;
	double **tmp_grid_curve_conc, **tmp_grid_curve;
	time_t start, end, start_true, end_true, time_true;
	node *new_node, *tmp;
	curve *tmp_real;
	tuple_true_neighbor *true_neigh;

	// Open input file
	FILE *query = fopen(query_file, "r");
  
  	if (query == NULL)
  	{
    	printf("Fopen: error opening query file!\n");
      	exit(1);
    }
	
	
	// Get diameter from query file
  	fgets(line, sizeof(line), query);
  	strtok(line,"\t");
  	
	data = strtok(NULL,"\n "); 
  	R = atof(data);
	
	// Allocate memory for array to store buckets indexes
	indexes = (int*) malloc(L * sizeof(int));
	
	if (indexes == NULL)
	{
		printf("Malloc: memory allocation error!\n");
		exit(3);
	}

  	while (fgets(line, sizeof(line), query) != NULL)
    {	
		curve *new_curve = (curve*) malloc(sizeof(curve));
		
		if (new_curve == NULL)
		{
			printf("Malloc: memory allocation error!\n");
			exit(3);
		}
		
		// Assign curve dimension
		new_curve->dimension = data_search->dimension;
		
    	data = strtok(line,"\t");
		
		// Assign curve's ID
		new_curve->ID_curve = atoi(data);
		
		stats[counter][0] = (double) (new_curve->ID_curve);

		data = strtok(NULL,"\t");
		
		// Assign the points of curve
		new_curve->points = atoi(data);
		
		// Allocate memory for array to store the points
		new_curve->coordinates = (double**) malloc((new_curve->points) * sizeof(double*));
		
		if (new_curve->coordinates == NULL)
		{
			printf("Malloc: memory allocation error!\n");
			exit(3);
		}

		// Parse every single point of curve
      	for (i = 0 ; i < new_curve->points ; i++)
        {	
			// Allocate memory for array to store each coordinate for every point (related to dimension {2,3,4})	
			new_curve->coordinates[i] = (double*) malloc(data_search->dimension * sizeof(double));
			
			if (new_curve->coordinates[i] == NULL)
			{
				printf("Malloc: memory allocation error!\n");
				exit(3);
			}
			
			// Cases for different dimensions
			switch(data_search->dimension)
			{
				case 2:
					data = strtok(NULL,",");
				
					// Assign 1st Coordinate
					data = data + 1;				
					new_curve->coordinates[i][0] = atof(data);
					
					// Assign 2st Coordinate
					data = strtok(NULL,")");				
					new_curve->coordinates[i][1] = atof(data);

					strtok(NULL," ");		
					break;
				case 3:
					data = strtok(NULL,",");			
				
					// Assign 1st Coordinate
					data = data + 1;
					new_curve->coordinates[i][0] = atof(data);

					// Assign 2st Coordinate
					data = strtok(NULL,",");
					new_curve->coordinates[i][1] = atof(data);
					
					// Assign 3st Coordinate
					data = strtok(NULL,")");
					new_curve->coordinates[i][2] = atof(data);				
					
					strtok(NULL," ");			
					break;
				case 4:
					data = strtok(NULL,",");
				
					// Assign 1st Coordinate
					data = data + 1;
					new_curve->coordinates[i][0] = atof(data);
				  
					// Assign 2st Coordinate
					data = strtok(NULL,",");
					new_curve->coordinates[i][1] = atof(data);
				  
					// Assign 3st Coordinate
					data = strtok(NULL,",");
					new_curve->coordinates[i][2] = atof(data);
					
					// Assign 4st Coordinate
					data = strtok(NULL,")");
					new_curve->coordinates[i][3] = atof(data);
					
					strtok(NULL," ");			
					break;
				default:
					printf("Wrond input dimension!\n");
					exit(1);	
			}	
        }

		// Temporary array to concatenate all grid curves
		tmp_grid_curve_conc = (double**) malloc((K * (new_curve->points)) * sizeof(double*));
		
		if (tmp_grid_curve_conc == NULL)
		{
			printf("Malloc: memory allocation error!\n");
			exit(3);
		}
		
		for (i = 0 ; i < (K * (new_curve->points)) ; i++)
		{
			tmp_grid_curve_conc[i] = (double*) malloc((new_curve->dimension) * sizeof(double));
			
			if (tmp_grid_curve_conc[i] == NULL)
			{
				printf("Malloc: memory allocation error!\n");
				exit(3);
			}
		}
		
		size_neigh_ID = 0;
		found_grid = 0;
		time_true = 0;
		min_distance = INFINITY;
		min_ID = -1;
		
		neigh_ID = (int*) malloc(size_neigh_ID * sizeof(int));

		start = time(NULL);
		
		// Fill L hashtables
		for (i = 0 ; i < L ; i++)
		{	
			
			start_true = time(NULL);
			
			// Find k grid curves
			for (j = 0 ; j < K ; j++)
			{			  
				// Find grid curve for each grid
				tmp_grid_curve = find_grid_curve(new_curve, data_search->d, data_search->array_grid_t[i][j]);
				
				for (z = (j * (new_curve->points)) ; z < ((j+1) * (new_curve->points)) ; z++)
					tmp_grid_curve_conc[z] = tmp_grid_curve[z - (j * (new_curve->points))];
			}
			
			// Create final grid vector	and store to struct _node_		
			new_node = create_hashtable_node(tmp_grid_curve_conc, new_curve, K);
			
			
			end_true = time(NULL);
			
			time_true += end_true - start_true;
			
			// Check which hash funtion to use (CLASSIC or LSH)
			if (hash_function == CLASSIC)
				index = hash_function_classic(new_node->grid_curve_1D, new_node->grid_points, 
											data_search->ht[i]->size, data_search->array_r_classic[i]);				
			else
				index = hash_function_LSH(new_node->grid_curve_1D, new_node->grid_points, data_search->ht[i]->size, 
											data_search->array_r_classic[i], i, data_search->array_LSH);

			// Fill the bucket's index for L's hashtable
			indexes[i] = index;
			
			tmp = data_search->ht[i]->table[index];
			
			// Search the whole list of this table[index]
			while (tmp != NULL)
			{
				// Only for curves with the same amount of grid points
				if (new_node->grid_points == tmp->grid_points)
				{	
					// Only for curves with same 1D grid curves
					if (compare_grid_curves(new_node->grid_curve_1D, tmp->grid_curve_1D, tmp->grid_points))
					{
						
						found_grid = 1;
						
						// Calculate distance
						if (metric_function == DFT)
							distance = frechet_distance(new_node->real_curve->coordinates, tmp->real_curve->coordinates, 
										new_node->real_curve->points, tmp->real_curve->points, tmp->real_curve->dimension);
						else
							distance = DTW_distance(new_node->real_curve->coordinates, tmp->real_curve->coordinates, 
										new_node->real_curve->points, tmp->real_curve->points, tmp->real_curve->dimension);	

										
						// Min LSH distance
						if (distance < stats[counter][1])
							stats[counter][1] = distance;
												
						// Max LSH distance
						if (distance > stats[counter][2])
							stats[counter][2] = distance;
						
						// Average LSH distance
						stats[counter][3] += distance;
						
						stats[counter][9]++;
						
						
						// R diameter neighbors						
						if (distance < R)
						{
							// Every single time that there is a "near-neighbor" ,memory is reallocated..
							neigh_ID = (int*) realloc(neigh_ID, (++size_neigh_ID) * sizeof(int));
							
							if (neigh_ID == NULL)
							{
								printf("Realloc: memory allocation error1!\n");
								exit(3);
							}
							
							//.. to store it
							neigh_ID[size_neigh_ID-1] = tmp->real_curve->ID_curve;
						}
						
					}
				}
				// Step to the next node in bucket
				tmp = (tmp->next);
			}	
		}	
		
		// If we dont find any same grid curve in buckets, check all the curves in these buckets
		if (!found_grid)
		{
			// For L hashtables
			for (i = 0 ; i < L ; i++)
			{
				// Assign for the i' hashtable the i' index of query curve
				tmp = data_search->ht[i]->table[indexes[i]];
				
				// For the whole bucket with this index
				while (tmp != NULL)
				{
					// Calculate distance
					if (metric_function == DFT)
						distance = frechet_distance(new_node->real_curve->coordinates, tmp->real_curve->coordinates, 
									new_node->real_curve->points, tmp->real_curve->points, tmp->real_curve->dimension);
					else
						distance = DTW_distance(new_node->real_curve->coordinates, tmp->real_curve->coordinates, 
									new_node->real_curve->points, tmp->real_curve->points, tmp->real_curve->dimension);
					
					
					// Nearest neighbor
					if (distance < min_distance)
					{
						min_distance = distance;
						min_ID = tmp->real_curve->ID_curve;
					}
					
					// Min LSH distance
					if (distance < stats[counter][1])
						stats[counter][1] = distance;
											
					// Max LSH distance
					if (distance > stats[counter][2])
						stats[counter][2] = distance;
					
					// Average LSH distance
					stats[counter][3] += distance;
					
					stats[counter][9]++;
					
											
					// R diameter neighbors						
					if (distance < R)
					{
						neigh_ID = (int*) realloc(neigh_ID, (++size_neigh_ID) * sizeof(int));
						
						if (neigh_ID == NULL)
						{
							printf("Realloc: memory allocation error2!\n");
							exit(3);
						}
						
						neigh_ID[size_neigh_ID-1] = tmp->real_curve->ID_curve;
					}
					
					// Step
					tmp = (tmp->next);
				}
			}
		}
			
		// Sort the array that contains R-neighbor's ID 
		quickSort(neigh_ID, 0, size_neigh_ID-1);
		
		// Remove duplicates from neigh_ID
		neigh_ID_no_duplicates = remove_duplicates(neigh_ID, size_neigh_ID);
		
		
		end = time(NULL);
		
		if ((end - start) < stats[counter][5])
			stats[counter][5] = (double) (end - start);
		
		if ((end - start) > stats[counter][6])
			stats[counter][6] = (double) (end - start);
		
		stats[counter][7] += (double) (end - start);
		

		start = time(NULL);
		
		// Find the true neighbor for the curve
		true_neigh = find_true_nearest_neighbor(data_search->ht[0], new_node, L, metric_function);
		
		end = time(NULL);
		
		
		// True distance
		stats[counter][4] = true_neigh->true_min_distance;
		
		// True time
		stats[counter][8] =  (end - start) + time_true;
				
		counter++;
		
		// Free new_node			
		tmp_real = new_node->real_curve;

		for (z = 0 ; z < tmp_real->points ; z++)
			free(tmp_real->coordinates[z]);
		free(tmp_real->coordinates);

		free(tmp_real);	
		
		free(new_node->grid_curve_1D);
		free(new_node);
		
		// Free neigh_ID_no_duplicates	
		free(neigh_ID_no_duplicates);
				
		// Free true_neigh
		free(true_neigh);
	}
	
	// Free indexes
	free(indexes);
	
	fclose(query);
	return;
}


// Function that checks if two grid curves are similar
int compare_grid_curves(double *grid_curve_1, double *grid_curve_2, int size)
{
	int i;
	
	for (i = 0 ; i < size ; i++)
		if(grid_curve_1[i] != grid_curve_2[i])
			return 0;
		
	return 1;
}


// Function that finds out the true nearest neighbor of the query curve
tuple_true_neighbor *find_true_nearest_neighbor(hashtable *ht, node *query_node, int L, int metric_function)
{
	int i, true_min_ID = -1;
	double distance, true_min_distance = INFINITY;
	node *tmp;
	

	for (i = 0 ; i < ht->size ; i++)
	{
		tmp = ht->table[i];

		while (tmp != NULL)
		{				
			// Calculate distance
			if (metric_function == DFT)
				distance = frechet_distance(query_node->real_curve->coordinates, tmp->real_curve->coordinates, 
							query_node->real_curve->points, tmp->real_curve->points, tmp->real_curve->dimension);
			else
				distance = DTW_distance(query_node->real_curve->coordinates, tmp->real_curve->coordinates, 
							query_node->real_curve->points, tmp->real_curve->points, tmp->real_curve->dimension);
				
				
			// Nearest true neighbor
			if (distance < true_min_distance)
			{
				true_min_distance = distance;
				true_min_ID = (tmp->real_curve->ID_curve);
			}

			// Step
			tmp = (tmp->next);
		}					
	}
	
	// Allocate memory for struct true_neigh_node that keeps true_min_ID, true_min_distance
	tuple_true_neighbor *true_neigh_node = (tuple_true_neighbor*) malloc(sizeof(tuple_true_neighbor));

	if (true_neigh_node == NULL)
	{
		printf("Malloc: memory allocation error!\n");
		exit(3);
	}

	// Fill in the fields
	true_neigh_node->true_min_ID = true_min_ID;
	true_neigh_node->true_min_distance = true_min_distance;

	return true_neigh_node;	
}


// Function that removes consecutive duplicates values from neigh_ID 
int *remove_duplicates(int *neigh_ID, int size_neigh_ID)
{
	int i, size_neigh_ID_no_dupl = 1;
	int *neigh_ID_no_duplicates;
	
	// Fill in the last position of array with -1 , as flag of the end of array
	neigh_ID = (int*) realloc(neigh_ID, (++size_neigh_ID) * sizeof(int));
		
	if (neigh_ID == NULL)
	{
		printf("Realloc: memory allocation error3!\n");
		exit(3);
	}
	
	neigh_ID[size_neigh_ID-1] = -1;


	// Allocate memory for new array that keeps our unique values ( no duplicates )
	neigh_ID_no_duplicates = (int*) malloc(size_neigh_ID_no_dupl * sizeof(int));
	
	if (neigh_ID_no_duplicates == NULL)
	{
		printf("Malloc: memory allocation error3!\n");
		exit(3);
	}
	
	neigh_ID_no_duplicates[0] = neigh_ID[0];
	
	// For the size of neigh_ID array .. 
	for (i = 1 ; i < size_neigh_ID ; i++)
	{
		// Only if there is a different neighbor's ID..
		if (neigh_ID_no_duplicates[size_neigh_ID_no_dupl-1] != neigh_ID[i])
		{
			//.. we realloc memory ..
			neigh_ID_no_duplicates = (int*) realloc(neigh_ID_no_duplicates, (++size_neigh_ID_no_dupl) * sizeof(int));
			
			if (neigh_ID_no_duplicates == NULL)
			{
				printf("Realloc: memory allocation error3!\n");
				exit(3);
			}
			
			//.. to store it
			neigh_ID_no_duplicates[size_neigh_ID_no_dupl-1] = neigh_ID[i];
		}
	}
		
	free(neigh_ID);
	
	return neigh_ID_no_duplicates;
}


// Function that returns the number of curves in query file
int number_of_curves_query(char *query_file)
{		
	int number_of_curves = 0;
	char line[LINE_LENGTH];
	char *data;
	
	FILE *query = fopen(query_file, "r");
	
	if (query == NULL)
  	{
    	printf("Fopen: error opening input file!\n");
      	exit(1);
    }
	
	// Skip first line
	fgets(line, sizeof(line), query);
	
	while(fgets(line, sizeof(line), query) != NULL)
		number_of_curves++;
	

	fclose(query);
	
	return number_of_curves;
}