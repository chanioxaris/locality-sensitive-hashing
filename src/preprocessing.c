#include "functions.h"

// Function that reads the data from the input file and creates the hashtables
database *preprocessing(char *input_file, int L, int K, int hash_function)
{
  	int i, j, z, dimension, index;
	int **array_r_classic;
	double d;
	double **array_grid_t, **tmp_grid_curve_conc, **tmp_grid_curve, **array_LSH;
	char line[LINE_LENGTH];
  	char *data;	
	node *new_node;		// Node for hashtable's linked list that contains all the info for each curve
	hashtable **ht;	
	
	dataset_info *dataset_information = (dataset_info*) malloc(sizeof(dataset_info));
	dataset_information = get_dataset_information(input_file);
	
	// Open input file
	FILE *input = fopen(input_file, "r");
	
	// Check if input file has opened successfully
  	if (input == NULL)
  	{
    	printf("Fopen: error opening input file!\n");
      	exit(1);
    }
	
	// Initialization of array that contains L hashtables
	ht = hashtable_init(L, dataset_information->number_of_curves);
	
	// Get dimension from input file
  	fgets(line, sizeof(line), input);
  	strtok(line,"\t"); 	
	data = strtok(NULL,"\n "); 
  	dimension = atoi(data);
	
	// Calculate step d for grid
	d = (double)(4 * dimension * (dataset_information->min_points))/DECIMAL_ACCURACY;
	
	// Generate array that contains transposition t for each grid
	array_grid_t = generate_array_grid_t(L, K, d);
		
	// Generate array that contains r values for linear combination
	array_r_classic = generate_array_r_classic(L, (dataset_information->max_points) * dimension);
	
	// Generate array that contains KVEC*L transpositions and the v vectors
	array_LSH = generate_array_LSH(L, (dataset_information->max_points) * dimension);

	// Free dataset_information
	free(dataset_information);
	
  	while (fgets(line, sizeof(line), input) != NULL)
    {	
		curve *new_curve = (curve*) malloc(sizeof(curve));
		
		if (new_curve == NULL)
		{
			printf("Malloc: memory allocation error!\n");
			exit(3);
		}
		
		// Assign curve dimension
		new_curve->dimension = dimension;
		
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
			new_curve->coordinates[i] = (double*) malloc(dimension * sizeof(double));
			
			if (new_curve->coordinates[i] == NULL)
			{
				printf("Malloc: memory allocation error!\n");
				exit(3);
			}
			
			// Cases for different dimensions
			switch(dimension)
			{
				case 2:
					data = strtok(NULL,",");
				
					// Assign 1st coordinate
					data = data + 1;				
					new_curve->coordinates[i][0] = atof(data);
					
					// Assign 2st coordinate
					data = strtok(NULL,")");				
					new_curve->coordinates[i][1] = atof(data);

					strtok(NULL," ");		
					break;
				case 3:
					data = strtok(NULL,",");			
				
					// Assign 1st coordinate
					data = data + 1;
					new_curve->coordinates[i][0] = atof(data);

					// Assign 2st coordinate
					data = strtok(NULL,",");
					new_curve->coordinates[i][1] = atof(data);
					
					// Assign 3st coordinate
					data = strtok(NULL,")");
					new_curve->coordinates[i][2] = atof(data);				
					
					strtok(NULL," ");			
					break;
				case 4:
					data = strtok(NULL,",");
				
					// Assign 1st coordinate
					data = data + 1;
					new_curve->coordinates[i][0] = atof(data);
				  
					// Assign 2st coordinate
					data = strtok(NULL,",");
					new_curve->coordinates[i][1] = atof(data);
				  
					// Assign 3st coordinate
					data = strtok(NULL,",");
					new_curve->coordinates[i][2] = atof(data);
					
					// Assign 4st coordinate
					data = strtok(NULL,")");
					new_curve->coordinates[i][3] = atof(data);
					
					strtok(NULL," ");			
					break;
				default:
					printf("Wrond input dimension!\n");
					exit(1);	
			}	
        }
	
		// Temporary array to concatenate K grid curves
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
		
		// Fill L hashtables
		for (i = 0 ; i < L ; i++)
		{
			// Find K grid curves
			for (j = 0 ; j < K ; j++)
			{			  
				// Find grid curve for every grid
				tmp_grid_curve = find_grid_curve(new_curve, d, array_grid_t[i][j]);
				
				// Fill tmp_grid_curve_conc array with each grid curve to accomplish concatenation
				for (z = (j * (new_curve->points)) ; z < ((j+1) * (new_curve->points)) ; z++)
					tmp_grid_curve_conc[z] = tmp_grid_curve[z - (j * (new_curve->points))];
			}
			
			// Create final grid vector	and store data to struct _node_		
			new_node = create_hashtable_node(tmp_grid_curve_conc, new_curve, K);
			
			// Check which hash funtion to use (CLASSIC or LSH)
			if (hash_function == CLASSIC)
				index = hash_function_classic(new_node->grid_curve_1D, new_node->grid_points, ht[i]->size, array_r_classic[i]);								
			else
				index = hash_function_LSH(new_node->grid_curve_1D, new_node->grid_points, ht[i]->size, array_r_classic[i], i, array_LSH);			
			
			// Insert new_node to hashtable 
			hashtable_insert(ht[i], index, new_node);
		}		
	}
	
	// Close input file
	fclose(input);
	
	// Allocate memory for struct that contains all needed info for the search operation
	database *data_search = (database*) malloc(sizeof(database));
	
	if (data_search == NULL)
	{
		printf("Malloc: memory allocation error!\n");
		exit(3);
	}
	
	// Fill data_search fields with data
	data_search->dimension = dimension;
	data_search->d = d;
	data_search->array_grid_t = array_grid_t;
	data_search->array_r_classic = array_r_classic;
	data_search->array_LSH = array_LSH;	
	data_search->ht = ht;
		
	return data_search;
}


// Function that returns min curve points, max curve points and number of curves in dataset
dataset_info *get_dataset_information(char *input_file)
{		
	int point, max_points = 0, min_points = INT_MAX, number_of_curves = 0;
	char line[LINE_LENGTH];
	char *data;
	
	dataset_info *info = (dataset_info*) malloc(sizeof(dataset_info));
	
	FILE *input = fopen(input_file, "r");
	
	if (input == NULL)
  	{
    	printf("Fopen: error opening input file!\n");
      	exit(1);
    }
	
	// Skip first line
	fgets(line, sizeof(line), input);
	
	while(fgets(line, sizeof(line), input) != NULL)
    {
		// Skip curve_ID
		strtok(line,"\t");
		
		// Get #points for this curve
		data = strtok(NULL,"\t");
		
		point = atoi(data);
		
		// Assign the new minimum points if data < old min_points
		if (point > max_points)
			max_points = point;		
		
		// Assign the new minimum points if data < old min_points
		if (point < min_points)
			min_points = point;
		
		number_of_curves++;
	}

	fclose(input);
	
	info->min_points = min_points;
	info->max_points = max_points;
	info->number_of_curves = number_of_curves;
	
	return info;
}


// Function that generates L*K 	transpositions for each grid 
double **generate_array_grid_t(int L, int K, double d)
{
	int i, j;
	double **array;
	
	// Allocate space for new array
	array = (double**) malloc(L * sizeof(double*));
	
	if (array == NULL)
	{
		printf("Malloc: memory allocation error!\n");
		exit(3);
	}
	
	for (i = 0 ; i < L ; i++)
	{
		array[i] = (double*) malloc(K * sizeof(double));
		
		if (array[i] == NULL)
		{
			printf("Malloc: memory allocation error!\n");
			exit(3);
		}
	}
	
	// Generate numbers that follow Gaussian distribution
	for (i = 0 ; i < L ; i++)
		for (j = 0 ; j < K ; j++)
			array[i][j] = d * rand_gaussian();
		
	return array;	
}


// Function that generates a random number following Gaussian Distribution
double rand_gaussian()
{
	double r, y1, y2;

	do
	{
		y1 = (rand() / (RAND_MAX + 1.0)) * 2.0;
		y2 = (rand() / (RAND_MAX + 1.0)) * 2.0;  
		r = y1*y1 + y2*y2;	
	}
	while (r >= 1);

	return y1 * pow(((pow(r,-2) - 1)/r) , 1/2);
}