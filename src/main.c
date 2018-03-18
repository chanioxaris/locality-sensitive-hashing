#include "functions.h"

int main(int argc, char* argv[])
{
  	// Initialize variables
  	int i, j, K = 2, L = 3, stats = 0, metric_function = 0, hash_function = 0, query_size;
	char input_file[PATH_LENGTH] = "NULL", query_file[PATH_LENGTH] = "NULL", output_file[PATH_LENGTH] = "NULL", user_input[16];
	double **stats_array;
	
	database *data;

    //Parse arguments from command line
    if (argc >= 1 || argc <= 16)									
    {
        for (i = 1 ; i <= (argc-1) ; i++)
        {
            if (!strcmp(argv[i], "-d")) // If flag is -d ..
            {
                stpcpy(input_file, argv[++i]); // Copy of argument's value to var
                continue;
            }
            if (!strcmp(argv[i], "-q")) // If flag is -q ..
            {
                stpcpy(query_file, argv[++i]); // Copy of argument's value to var
                continue;
            }
            if (!strcmp(argv[i], "-k"))	// If flag is -k ..
            {
                K = atoi(argv[++i]);	// Convertion from char to int and var assignment
                continue;
            }
            if (!strcmp(argv[i], "-L")) // If flag is -k ..
            {
                L = atoi(argv[++i]);	// Convertion from char to int and var assignment
                continue;
            }
            if (!strcmp(argv[i], "-o")) // If flag is -o ..
            {
                stpcpy(output_file, argv[++i]); // Copy of argument's value to var
                continue;
            }
            if (!strcmp(argv[i], "-stats")) // If flag is -stats ..
            {
                stats = 1;	// Set this var True
                continue;
            }
            if (!strcmp(argv[i], "-function")) // If flag is -function ..
            {
                if (!strcmp(argv[++i], "DFT")) // Check if metric is "DFT"
                {
                    metric_function = DFT; 
                }
                else if (!strcmp(argv[i], "DTW")) // or if it's "DTW"
                {
                    metric_function = DTW; 
                }
                else	//Otherwise .. the input is wrong!
                {
                    printf("Wrong metric function argument! \n");
                    return -1;
                }

                continue;
            }
            if (!strcmp(argv[i], "-hash")) // If flag is -hash ..
            {
                if (!strcmp(argv[++i], "classic")) // Check if hash function is "classic"
                {
                    hash_function = CLASSIC;  
                }
                else if (!strcmp(argv[i], "probabilistic")) // Check if next argument is "probabilistic"
                {
                    hash_function = PROB; 
                }
                else	//Otherwise .. the input is wrong!
                {
                    printf("Wrong hash function argument! \n");
                    return -2;
                }

                continue;
            }	

            //if flag of any argument is wrong .. exit!
            printf("Wrong input arguments! \n");
            return -3;
            }
        }
    else	//if number of arguments is wrong .. exit!
    {
        printf("Wrong number of arguments! \n");
        return -4;	
    }


    //Parse input from user
    if (!strcmp(input_file, "NULL"))
    {
        printf("Please insert the path to input file: ");
        scanf("%s", input_file);
    } 

    if (!metric_function)
    {
        printf("Please insert the metric function: ");
        scanf("%s", user_input);

        if (!strcmp(user_input, "DFT")) // Check if metric is "DFT"
        {
            metric_function = DFT; 
        }
        else if (!strcmp(user_input, "DTW")) // or if it's "DTW"
        {
            metric_function = DTW; 
        }
        else	//Otherwise .. the input is wrong!
        {
            printf("Wrong metric function argument! \n");
            return -1;
        }      
    } 

    if (!hash_function)
    {
        printf("Please insert the hash function: ");
        scanf("%s", user_input);

        if (!strcmp(user_input, "classic")) // Check if hash function is "classic"
        {
            hash_function = CLASSIC;  
        }
        else if (!strcmp(user_input, "probabilistic")) // Check if next argument is "probabilistic"
        {
            hash_function = PROB; 
        }
        else	//Otherwise .. the input is wrong!
        {
            printf("Wrong hash function argument! \n");
            return -2;
        }
    } 
		
	// Seed random number generator 
	srand(time(NULL));

	if (!stats)
	{
		// Preprocessing
		printf("Creating search structure...\n");
		data = preprocessing(input_file, L, K, hash_function);

		// Iterate until user exits the program  
		do
		{
			if (!strcmp(query_file, "NULL"))
			{
				// Εnter the path of the query file
				printf("Please insert the path to query file: ");
				scanf("%s", query_file);
			} 
			
			if (!strcmp(output_file, "NULL"))
			{
				// Εnter the path of the output file
				printf("Please insert the path to output file: ");
				scanf("%s", output_file);
			} 
			   
			// Search
			printf("Entering search...\n");
			search(query_file, output_file, data, L, K, hash_function, metric_function);


			// Ask user for another search 
			stpcpy(query_file, "NULL");
			stpcpy(output_file, "NULL");

			printf("Do you want to proceed with another search? (Y/N): ");
			scanf("%s", user_input);
		 
		}while(!strcmp(user_input, "Y")); 
		
		
		// Free data struct
		for (i = 0 ; i < L ; i++)
			free(data->array_grid_t[i]);
		free(data->array_grid_t);

		for (i = 0 ; i < (KVEC*L) ; i++)
			free(data->array_LSH[i]);
		free(data->array_LSH);
		
		for (i = 0 ; i < L ; i++)
			free(data->array_r_classic[i]);
		free(data->array_r_classic);
		
		hashtable_destroy(data->ht, L);
		
		free(data);
	}
	else
	{
		// Iterate until user exits the program  
		do
		{
			if (!strcmp(output_file, "NULL"))
			{
				// Εnter the path of the output file
				printf("Please insert the path to output file: ");
				scanf("%s", output_file);
			} 

			if (!strcmp(query_file, "NULL"))
			{
				// Εnter the path of the query file
				printf("Please insert the path to query file: ");
				scanf("%s", query_file);
			} 
			
			query_size = number_of_curves_query(query_file);
			
			stats_array = (double**) malloc(query_size * sizeof(double*));
			
			if (stats_array == NULL)
			{
				printf("Malloc: memory allocation error!\n");
				exit(3);
			}
			
			for (i = 0 ; i < query_size ; i++) 
			{
				stats_array[i] = (double*) malloc(10 * sizeof(double));
				
				if (stats_array[i] == NULL)
				{
					printf("Malloc: memory allocation error!\n");
					exit(3);
				}
			}
							
			// Initialize array with stats
			for (i = 0 ; i < query_size ; i++) 
			{
				stats_array[i][1] = INFINITY;
				stats_array[i][2] = -INFINITY;
				stats_array[i][3] = 0.0;
				stats_array[i][4] = INFINITY;
				stats_array[i][5] = INFINITY;
				stats_array[i][6] = -INFINITY;
				stats_array[i][7] = 0.0;
				stats_array[i][8] = INFINITY;
				stats_array[i][9] = 0.0;
			}
				
			// Iterate 100 times to gain stats
			for (i = 0 ; i < ITERATIONS ; i++)
			{
				// Preprocessing
				printf("Creating search structure...\n");
				data = preprocessing(input_file, L, K, hash_function);
				 
				// Search
				printf("Entering search...\n");
				search_stats(query_file, data, stats_array, L, K, hash_function, metric_function);	

				
				// Free data struct
				for (j = 0 ; j < L ; j++)
					free(data->array_grid_t[j]);
				free(data->array_grid_t);

				for (j = 0 ; j < (KVEC*L) ; j++)
					free(data->array_LSH[j]);
				free(data->array_LSH);

				for (j = 0 ; j < L ; j++)
					free(data->array_r_classic[j]);
				free(data->array_r_classic);

				hashtable_destroy(data->ht, L);

				free(data);
			}
			
			output_stats(output_file, stats_array, metric_function, hash_function, query_size);
			
			// Free memory of stats_array			
			for (i = 0 ; i < query_size ; i++) 
				free(stats_array[i]);			
			free(stats_array);
			
					  
			// Ask user for another search 
			stpcpy(query_file, "NULL");	
			stpcpy(output_file, "NULL");

			printf("Do you want to proceed with another search? (Y/N): ");
			scanf("%s", user_input);
					  					  
		}while(!strcmp(user_input, "Y")); 				  
	}
      
	return 0; 
}