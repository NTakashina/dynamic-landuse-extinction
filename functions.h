//
//  functions.h
//

#ifndef _functions_h
#define _functions_h

// Species list structure
typedef struct range_list{
    int species_id;          // Species ID
    int sar;
    int ear;
    int ear2;                // EAR count for 95% habitat loss
    int ear3;                // EAR count for 90% habitat loss
    int area;                // Range size
    int sp_map[101][101];
    
    int species_id_backup;   // Species ID
    int area_backup;         // Range size
    int sp_map_backup[101][101];
    
    struct range_list *next;
}species;

species *initial_list(species *, int , int imax, int, int [][imax+1]);

void list_recovery_initial_condition(species * , int);
void list_input(species *, double , int , double , int);
void list_sar_ear(species *, int , double , double *[] , int *, int *);
int pa_ripley_lr(int , double , int *[], double , int , double *);
int pa_paircorr_g(int , double , int *[], double , int , double *);
void free_list(species *);
    
// 2D PDE solver
double rk_2d( int , int , int , int , double , double , double , double , double , double *[], double **[] );
void PDE(species *, int , int , double , double , double , double , double , double , double , double *[], int *[], double **[], int , int , int , int **[], double *[], int *[], int , int , double , double );

// Geographic range and converted area
void urban_pa_area(int , double *[], int *[], double *, double * );

// Abundance
double abundance_2d( double *[], int , double );

// Random number generators
int rand_po(double);
double rand_gauss(double , double);

// Data
void data_connection(int **[], double *[], int *[], int , int, double , int , double , double );
void ripley_connection(double *[], int [], int , int , double );
void paircorr_connection(double *[], int [], int , int , double );

// Memory
void release_memory_int_2d(int *[], int);
void release_memory_int_3d(int **[], int , int);
void release_memory_double_2d(double *[], int);
void release_memory_double_3d(double **[], int , int);

#endif
