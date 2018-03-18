#include "functions.h"

// Function that appends to output file if no -stats was requested
void output_nonstats(char *output_file, int curve_ID, int metric_function, int hash_function, int found_grid,
                     int ID_LSH, int ID_true, double dist_LSH, double dist_true, int *neighbors)
{
  	int i;
  
  	// Open the output_file to append new information
	FILE *output = fopen(output_file, "a");
  	
  	// Check if the file opened successfully
  	if (output == NULL)
  	{
		printf("Fopen: error opening output file!\n");
      	exit(1);
    }
   
  	// Let's start writing in file..
	fprintf(output, "Query: %d\n", curve_ID);
	
  	if (metric_function == DFT)
    	fprintf(output, "DistanceFunction: DFT\n");
    else
        fprintf(output, "DistanceFunction: DTW\n");
   
  	if (hash_function == CLASSIC)
    	fprintf(output, "HashFunction: Classic\n");
    else
        fprintf(output, "HashFunction: Probabilistic\n");

  	if (found_grid)
      	fprintf(output, "FoundGridCurve : True\n");
	else
      	fprintf(output, "FoundGridCurve : False\n");

    fprintf(output, "LSH Nearest Neighbor: %d\n", ID_LSH);
  
    fprintf(output, "True Nearest Neighbor: %d\n", ID_true);
  
  	fprintf(output, "distanceLSH: %.16f\n", dist_LSH);
  
    fprintf(output, "distanceTrue: %.16f\n", dist_true);
  
    fprintf(output, "R-nearest neighbors:\n");  
  	for (i = 0 ; neighbors[i] != -1 ; i++)
      	fprintf(output, "%d\n", neighbors[i]);
    
  	fprintf(output, "\n\n\n"); 
  
  	// When the writing is finished, we close the file!
	fclose(output);
  
	return;
}


// Function that appends to output file if -stats was requested
void output_stats(char *output_file, double **stats, int metric_function, int hash_function, int query_size)
{
  	int i;
  
  	// Open the output_file to append new information
	FILE *output = fopen(output_file, "a");
  	
  	// Check if the file opened successfully
  	if (output == NULL)
  	{
    	printf("Fopen: error opening output file!\n");
      	exit(1);
    }
    
	for (i = 0 ; i < query_size ; i++)
	{
		// Let's start writing in file..
		fprintf(output, "Query: %d\n", (int) stats[i][0]);
		
		if (metric_function == DFT)
			fprintf(output, "DistanceFunction: DFT\n");
		else
			fprintf(output, "DistanceFunction: DTW\n");
	   
		if (hash_function == CLASSIC)
			fprintf(output, "HashFunction: Classic\n");
		else
			fprintf(output, "HashFunction: Probabilistic\n");
		

		fprintf(output, "|minDistanceLSH - distanceTrue|: %f\n", fabs(stats[i][1] - stats[i][4]));
		fprintf(output, "|maxDistanceLSH - distanceTrue|: %f\n", fabs(stats[i][2] - stats[i][4])); 
		fprintf(output, "|avgDistanceLSH - distanceTrue|: %f\n", fabs(((double) stats[i][3]/(double) stats[i][9]) - stats[i][4]));
	 
		fprintf(output, "tLSHmin: %d seconds\n", (int) stats[i][5]);
		fprintf(output, "tLSHmax: %d seconds\n", (int) stats[i][6]);
		fprintf(output, "tLSHavg: %.2f seconds\n", (double) stats[i][7]/ITERATIONS);
		fprintf(output, "tTrue: %d seconds\n", (int) stats[i][8]);
	  
		fprintf(output, "\n\n\n");
	}
  
  	// When the writing is finished, we close the file!
	fclose(output);
  
	return;
}