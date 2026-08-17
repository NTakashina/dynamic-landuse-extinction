//
//  urban_multi.c
//
//  Baseline implementation of the spatially explicit land-use model.
//
//  Land-use dynamics are solved on a 2D grid using a fourth-order
//  Runge-Kutta (RK4) time integration scheme.
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "functions.h"
#include "mt64.h"

#define XGRID 100
#define XMAX 1.0          // Maximum coordinate
#define XMIN 0.0          // Minimum coordinate
#define STEP_X 100        // Step size on x-axis
#define TMAX 100000
#define DT 0.005          // Step size of time
#define P_CYCLE 5         // Frequency to plot
#define SP 1000.0
#define PRINT 0           // 0 -> DO NOT PRINT DISTRIBUTION; 1 -> PRINT
#define TRIAL 100         // # of trials
#define SCENARIO 9        // # of scenarios
#define RIPLEY_DR 0.02
#define RIPLEY_R_MAX 0.45
#define PAIRCORR_DR 0.03
#define PAIRCORR_R_MAX 0.70

int main(void)
{

    FILE *fq;
    fq = fopen("data_ear.dat","w"), fclose(fq);
    
    // MT initialization
    init_genrand64((unsigned) time(NULL));
    
    species dummy_head;              // Dummy head node; stores no species data
    species *ptr;                    // Working pointer
    
    ptr = &dummy_head;
    ptr->next = NULL;
    
    int ***final_data;
    int  **final_cover;
    int  *ripley_count;
    int  *paircorr_count;
    double **u_x;        // Land-use map: u_x[i][j] > 0 indicates transformed area
    int    **u_p;        // Protected area map: u_p[i][j] = 1 indicates protected area
    double **data_time;
    double **ripley_sum;
    double **paircorr_sum;
    double ***u_matrix;  // Runge-Kutta stage values for the land-use field
    
    int i, i2, j, step = STEP_X, ripley_bin = (int)floor(RIPLEY_R_MAX / RIPLEY_DR + 1.0e-9), paircorr_bin = (int)floor(PAIRCORR_R_MAX / PAIRCORR_DR + 1.0e-9);
    int ripley_valid, paircorr_valid;
    
    // Biological parameters
    double r_x = 1.0;          // Growth rate
    double K = 10.0;           // Carrying capacity
    double Dx = 0.005;         // Diffusion coefficient
    double Dy = 0.005;         // Diffusion coefficient
    int dist_scenario = 0;     // Species-center distribution: 0 -> Uniform; 1 -> Normal
    
    // Land-use parameter
    int step_urban;             // Time-step interval between land-use events
    int step_pa = 1;            // Time-step interval for adding protected cells
    double lu_level_pa_est;     // Land-use fraction at which PA establishment starts
    double max_pa_frac = 0.05;  // Maximum fraction of protected area
    
    // NUMERICAL SETTINGS
    double dw = (XMAX - XMIN) / (double)XGRID;   // Grid spacing for the 2D PDE simulation
    double dw2 = (XMAX - XMIN) / (double)STEP_X; // Bin width for EAR/SAR output
    
    printf("Dx * DT / (dw*dw) = %f\n", Dx * DT / (dw * dw));
    
    // MEMORY ALLOCATION
    final_data = (int ***)calloc(SCENARIO, sizeof(int **)), final_cover = (int **)calloc(SCENARIO, sizeof(int *)), data_time = (double **)calloc(SCENARIO, sizeof(double *));
    ripley_sum = (double **)calloc(SCENARIO, sizeof(double *)), ripley_count = (int *)calloc(SCENARIO, sizeof(int));
    paircorr_sum = (double **)calloc(SCENARIO, sizeof(double *)), paircorr_count = (int *)calloc(SCENARIO, sizeof(int));
    for (i = 0; i < SCENARIO; i++) {
        final_data[i] = (int **)calloc(STEP_X, sizeof(int *));
        final_cover[i] = (int *)calloc(STEP_X, sizeof(int));
        data_time[i] = (double *)calloc(20 + 1, sizeof(double));
        ripley_sum[i] = (double *)calloc(ripley_bin, sizeof(double));
        paircorr_sum[i] = (double *)calloc(paircorr_bin, sizeof(double));
        for(j = 0; j < STEP_X; j++){
            final_data[i][j] = (int *)calloc(4, sizeof(int));
        }
    }
    
    u_matrix = (double ***)calloc(XGRID + 1, sizeof(double **));
    for(i = 0; i < XGRID + 1; i++) {
        u_matrix[i] = (double **)calloc(XGRID + 1, sizeof(double *));
        for (j = 0; j < XGRID + 1; j++) {
            u_matrix[i][j] = (double *)calloc(5, sizeof(double));
        }
    }
        
    // Simulation starts here
    for(i2 = 0; i2 < TRIAL; i2++){
        printf("i2 = %d PROGRESS = %f%%\n", i2, 100.0 * (double)i2 / (double)TRIAL);
    
        // Initial condition
        list_input(&dummy_head, SP, XGRID, dw, dist_scenario);
        //histo_list_area(&dummy_head,XGRID,dw);  // Generate a histogram of range sizes (optional)
                       
        for(i = 0; i < SCENARIO; i++){
            // MEMORY ALLOCATION
            printf("i = %d\n", i);
            u_x = (double **)calloc(XGRID + 1, sizeof(double *));
            u_p = (int **)calloc(XGRID + 1, sizeof(int *));
            for (j = 0; j < XGRID + 1; j++) {
                u_x[j] = (double *)calloc(XGRID + 1, sizeof(double));
                u_p[j] = (int *)calloc(XGRID + 1, sizeof(int));
            }
                        
            if(i == 0) {step_urban = 999999999, lu_level_pa_est = 0.2;}
            else if (i == 1) {step_urban = 100, lu_level_pa_est = 0.2;}
            else if (i == 2) {step_urban = 10, lu_level_pa_est = 0.2;}
            else if (i == 3) {step_urban = 999999999, lu_level_pa_est = 0.4;}
            else if (i == 4) {step_urban = 100, lu_level_pa_est = 0.4;}
            else if (i == 5) {step_urban = 10, lu_level_pa_est = 0.4;}
            else if (i == 6) {step_urban = 999999999, lu_level_pa_est = 1.0;}  // No PA
            else if (i == 7) {step_urban = 100, lu_level_pa_est = 1.0;}
            else if (i == 8) {step_urban = 10, lu_level_pa_est = 1.0;}
            
            PDE(&dummy_head,TMAX,XGRID,DT,Dx,Dy,r_x,K,dw,dw2,u_x,u_p,u_matrix,PRINT,step_urban,step_pa,final_data,data_time,final_cover,i,step,lu_level_pa_est,max_pa_frac);
            {
                double ripley_lr[ripley_bin];
                double paircorr_g[paircorr_bin];
                ripley_valid = pa_ripley_lr(XGRID, dw, u_p, RIPLEY_DR, ripley_bin, ripley_lr);
                paircorr_valid = pa_paircorr_g(XGRID, dw, u_p, PAIRCORR_DR, paircorr_bin, paircorr_g);
                if(ripley_valid == 1 && paircorr_valid == 1){
                    ripley_count[i]++;
                    for(j = 0; j < ripley_bin; j++) ripley_sum[i][j] += ripley_lr[j];
                    paircorr_count[i]++;
                    for(j = 0; j < paircorr_bin; j++) paircorr_sum[i][j] += paircorr_g[j];
                }
            }
            list_recovery_initial_condition(&dummy_head, XGRID);  // Reset species ranges to the initial condition
            //list_print(&dummy_head,XGRID,dw);
            
            release_memory_double_2d(u_x,XGRID + 1);
            release_memory_int_2d(u_p,XGRID + 1);
        }
    }
    
    // Final data input
    data_connection(final_data,data_time,final_cover,SCENARIO,STEP_X,dw2,TRIAL,max_pa_frac,lu_level_pa_est);
    ripley_connection(ripley_sum, ripley_count, SCENARIO, ripley_bin, RIPLEY_DR);
    paircorr_connection(paircorr_sum, paircorr_count, SCENARIO, paircorr_bin, PAIRCORR_DR);
    
    // Memory release
    release_memory_int_3d(final_data, SCENARIO, STEP_X);
    release_memory_int_2d(final_cover, SCENARIO);
    release_memory_double_2d(data_time, SCENARIO);
    release_memory_double_2d(ripley_sum, SCENARIO);
    release_memory_double_2d(paircorr_sum, SCENARIO);
    free(ripley_count);
    free(paircorr_count);
    release_memory_double_3d(u_matrix, XGRID + 1, XGRID + 1);
    free_list(&dummy_head);
    
    return 0;
    
}
