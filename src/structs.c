typedef struct _curve_ {
	int ID_curve;
	int dimension;
	int points;					// Number of curve points
	double **coordinates;
}curve;


typedef struct _node_ {
	curve *real_curve;
	int grid_points;
	double *grid_curve_1D;
	struct _node_ *next;
}node;


typedef struct _hashtable_ {
	int size;
	struct _node_ **table;
}hashtable;


typedef struct _dataset_info_ {
	int min_points;
	int max_points;
	int number_of_curves;
}dataset_info;


typedef struct _database_ {
	int dimension;
	double d;
	double **array_grid_t;
	double **array_LSH;
	int **array_r_classic;
	struct _hashtable_ **ht;
}database;


typedef struct _tuple_true_neighbor_ {
	int true_min_ID;
	double true_min_distance;
}tuple_true_neighbor;
