#include "functions.h"
#include "hashtable.h"
#include "preprocessing.h"

// Function that initialize an array of L hashtables
hashtable **hashtable_init(int L, int dataset_size) {
    int i, j;

	// Allocate memory for array that contains L pointers to struct _hashtable_
    hashtable **hashtable_array = (hashtable**) malloc(L * sizeof(hashtable*));

	if (hashtable_array == NULL) {
		printf("Malloc: memory allocation error!\n");
		exit(3);
	}

	// For L hashtables
    for (i = 0 ; i < L ; i++) {
		// Allocate memory for hashtable struct
        hashtable_array[i] = (hashtable*) malloc(sizeof(hashtable));

		if (hashtable_array[i] == NULL) {
			printf("Malloc: memory allocation error!\n");
			exit(3);
		}

		// Set the size of hashtable
        hashtable_array[i]->size = dataset_size/4;

		// Allocate memory for hashtable's tables
        hashtable_array[i]->table = (node**) malloc((hashtable_array[i]->size) * sizeof(node*));

		if (hashtable_array[i]->table == NULL) {
			printf("Malloc: memory allocation error!\n");
			exit(3);
		}

		// Initialize each table to NULL
        for (j = 0 ; j < (hashtable_array[i]->size) ; j++)
            hashtable_array[i]->table[j] = NULL;
    }

    return hashtable_array;
}


// Function that inserts a new node in a hashtable
void hashtable_insert(hashtable *ht, int index, node *new_node) {
    if (ht->table[index] == NULL)
        ht->table[index] = new_node;
    else {
        node *tmp = ht->table[index];
		new_node->next = tmp;
        ht->table[index] = new_node;
    }

    return;
}


// Free memory for every field of hashtable
void hashtable_destroy(hashtable **ht, int L) {
	int i, j, z;

	node *tmp_node;
	curve *tmp_real;

	for (i = 0 ; i < L ; i++) {
		for (j = 0 ; j < ht[i]->size ; j++) {
			tmp_node = ht[i]->table[j];

			while (tmp_node != NULL) {
				// Free once real curve struct
				if (i == 0) {
					tmp_real = tmp_node->real_curve;

					for (z = 0 ; z < tmp_real->points ; z++)
						free(tmp_real->coordinates[z]);
					free(tmp_real->coordinates);

					free(tmp_real);
				}

				free(tmp_node->grid_curve_1D);

				tmp_node = tmp_node->next;
			}
		}
	}

	for (i = 0 ; i < L ; i++)
		free(ht[i]);
	free(ht);

	return;
}


/* Function that generates the r for the classic hash function
		- size: max points of real curves * dimension */
int **generate_array_r_classic(int L, int size) {
    int i, j;
    int **array_r_classic;

	// Allocate memory for the array that keeps L arrays of random r variables
	array_r_classic = (int**) malloc(L * sizeof(int*));

	if (array_r_classic == NULL) {
		printf("Malloc: memory allocation error!\n");
		exit(3);
	}

	// For L hashtables
    for (i = 0 ; i < L ; i++) {
		// Allocate memory for each array of random r variables
		array_r_classic[i] = (int*) malloc(size * sizeof(int));

		if (array_r_classic[i] == NULL) {
			printf("Malloc: memory allocation error!\n");
			exit(3);
		}

		// Fill in the arrays with random numbers in range [0,999]
		for (j = 0 ; j < size ; j++)
			array_r_classic[i][j] = rand() % 1000;
	}

    return array_r_classic;
}


// Function for classic hashing
int hash_function_classic(double *grid_curve_1D, int grid_points, int ht_size, int *array_r_classic) {
	int i, sum = 0;
	unsigned int M = UINT_MAX - 4;   //= 2^32 - 5

	// sums = sum + (ri * pi) mod M , i = 1,2,3...grid_points
	for(i = 0 ; i < grid_points ; i++) {
    	sum += (int) (array_r_classic[i] * grid_curve_1D[i]) % M;
      	sum = sum % M;
    }

	// Return the absolute value of sum mod TableSize , because there
	// possibility to return negative index
    return abs(sum % ht_size);
}


// Function that generates the v for the LSH hash function
double **generate_array_LSH(int L, int size) {
    int i, j;
    double **array_LSH;

	// Allocate memory for the array that keeps KVEC*L arrays of random v,t variables
	// for this project 3 * 3 array that keeps v vectors and t variables (one for each ht)
	array_LSH = (double**) malloc((KVEC*L) * sizeof(double*));

	if ((array_LSH) == NULL) {
		printf("Malloc: memory allocation error!\n");
		exit(3);
	}

	// For (KVEC*L) ..
    for (i = 0 ; i < (KVEC*L) ; i++) {
		// Allocate memory for the first value of our array ( that will be t variable )
		array_LSH[i] = (double*) malloc((size + 1) * sizeof(double));

		if (array_LSH[i] == NULL) {
			printf("Malloc: memory allocation error!\n");
			exit(3);
		}

		// Fill in the value of t ( uniformly chosen )
		array_LSH[i][0] = (double) rand_gaussian() * W;	// Generate LSH's t in range [0, W)

		// Generate uniformly distributed vector v for LSH hashing
		for (j = 1 ; j < (size + 1) ; j++)
			array_LSH[i][j] = rand_gaussian();
	}

    return array_LSH;
}


/* Function for LSH hashing
		- ht_number: number of hashtable {1,2,...,L} */
int hash_function_LSH(double *grid_curve_1D, int grid_points, int ht_size, int *array_r_classic, int ht_number, double **array_LSH) {
	int i, j, sum;
	double sum_kvec;
	unsigned int M = UINT_MAX - 4;
	int *vector_LSH;

	// Allocate memory for array to store new vector after LSH hashing, kvec-dimension
	// vector_LSH is [(p*v1 + t1)/w, (p*v2 + t2)/w, (p*v3 + t3)/w]
	vector_LSH = (int*) malloc(KVEC * sizeof(int));

	if (vector_LSH == NULL) {
		printf("Malloc: memory allocation error!\n");
		exit(3);
	}

	// LSH hashing
	for(i = 0; i < KVEC; i++) {
		sum_kvec = 0.0;

		for(j = 0; j < grid_points; j++)
			// Thats the linear combination of vi * pi , i = 1,2,...,grid_points
			sum_kvec += (array_LSH[(ht_number*KVEC)+i][j+1] * grid_curve_1D[j]);

		// Store (sum + t) / w value , for each hashtable ( of L hashtables) to each vector_LSH variable
		vector_LSH[i] = (int) floor((sum_kvec + array_LSH[(ht_number*KVEC)+i][0])/W);
	}

	// For KVEC ..
	for(i = 0; i < KVEC; i++) {
    	sum += (int) (array_r_classic[i] * vector_LSH[i]) % M;
      	sum = sum % M;
    }

	free(vector_LSH);

    return abs(sum % ht_size);
}


// Auxiliary function that prints the content of a hashtable
void hashtable_print(hashtable *ht) {
	int i,j;
	node *tmp;

	for (i = 0 ; i < ht->size ; i++) {
		printf("INDEX: %d\n", i);
		tmp = ht->table[i];

		while (tmp != NULL) {
			printf("%d %d %d %d\n",tmp->real_curve->ID_curve, tmp->real_curve->points, tmp->real_curve->dimension, tmp->grid_points);

            for (j = 0 ; j < tmp->grid_points ; j++)
				printf("%.16f ", tmp->grid_curve_1D[j]);
			printf("\n");

			tmp = (tmp->next);
		}
	}
    
	return;
}
