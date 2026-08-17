//
//  functions.c
//
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include "functions.h"
#include "mt64.h"
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#define EPS 1.0e-2
#define EPS2 1.0e-12
#define EPS3 0.01 // I: threshold of land-use intensity
#define BIN_MAX 20
#define CONV_SAMPLE_TIME 0.5
#define CONV_WINDOW_TIME 10.0
#define CONV_MIN_TIME 20.0
#define CONV_TOL_AREA 2.0e-3
#define CONV_TOL_POP 2.0e-3
#define CONV_WINDOW_COUNT 21
#define COVERAGE_LEVEL 0.8

/*
 Operations of structure
 */
species *list_initial(species *ptr2, int sp_id, int imax, int sp_area, int range_s[][imax+1])
{
    int i, j;
    
    species *q;        // Pointer to the newly added node
    if((q = (species *) malloc( sizeof(species))) == NULL) {
        printf("malloc error\n");
        exit(EXIT_FAILURE);
    }
        
    // Initialize structure members
    q->species_id = sp_id;
    q->sar = 0;
    q->ear = 0;
    q->ear2 = 0; // 5%
    q->ear3 = 0; // 10%
    q->area = sp_area;
    
    // Copy the current state to the backup fields
    q->species_id_backup = sp_id;
    q->area_backup = sp_area;
    
    for(i = 0; i <= imax; i++) {
        for(j = 0; j <= imax; j++) {
            q->sp_map[i][j] = range_s[i][j];
            q->sp_map_backup[i][j] = range_s[i][j];
        }
    }
    
    //printf("id=%d area=%d area(95%%)=%d area(90%%)=%d\n", q->species_id, q->area_backup, (int)((double)(q->area_backup) * 0.05), (int)((double)(q->area_backup) * 0.1));
    
    // Link the new node to the end of the list
    ptr2->next = q;           // Redirect the last node to the new node
    q->next = NULL;           // Mark the new node as the end of the list
    ptr2 = q;
    
    return ptr2;              // Return the last node
    
}

void list_recovery_initial_condition(species *ptr2, int imax)
{
    
    int i, j;
    
    species *q = ptr2;        // Set pointer to list head
    q = q->next;
        
    while ( q != NULL ) {
        // Restore the initial values
        q->species_id = q->species_id_backup;
        q->sar = 0;
        q->ear = 0;
        q->ear2 = 0;
        q->ear3 = 0;
        q->area = q->area_backup;
        
        for( i = 0; i <= imax; i++ ) {
            for( j = 0; j <= imax; j++ ) {
                q->sp_map[i][j] = q->sp_map_backup[i][j];
            }
        }
        q = q->next;
    }
}

void list_input(species *ptr2, double SP, int imax, double dx, int dist_scenario){ // Define species distributions
    
    //FILE *fq;
    //fq = fopen("data_sp_center.dat","w");
    
    const gsl_rng_type * T; // Random number generator type
    gsl_rng * r;

    // Choose a generator from GSL_RNG_TYPE
    gsl_rng_env_setup();
    T = gsl_rng_default;  // default is mt19937
    r = gsl_rng_alloc(T); // Generate instance
    
    gsl_rng_set(r,time(NULL)); // Initialize the generator
    
    // range_s[i][j]: 0 -> absence; 1 -> presence
    
    int i, j, k, ix, iy;
    int status;                // Check whether the species center overlaps the domain
    int tot_sp;                // Total species in the whole area
    int count_sp = 0;          // Count species intersecting [0,1] x [0,1]; add them to the list
    int count_area;            // Count the geographic range area of a species
    
    int range_s[imax+1][imax+1];
    
    double dist, dist_x0, dist_x1, dist_y0, dist_y1, rd, x, y, xc, yc, p, variance_gauss;
    
    // Define the radius of the geographic range, rd
    // Beta distribution
    double a = 2.0;
    double b = 10.0;
    double ave_rd = a / (a + b);
    
    int outer_i = 0.5 + 2.0 * ave_rd / dx; // Outer extension of [0,1]*[0,1] that adds an extra 2*ave_rd grids to each side
    if(outer_i > imax) outer_i = imax;     // Set maximum outer_i due to the range of a beta dist. [0,1]
        
    if(SP < 250.0){tot_sp = rand_po(SP);}            // Number of species ~ Poisson(mu_s)
    else {tot_sp = (int)(rand_gauss(SP, SP) + 0.5);} // Number of species ~ N(mu_s, mu_s); Gaussian approximation
        
    if(tot_sp != 0)
    for(k = 1; k <= tot_sp; k++){ // SPECIES;
        
        status = 0;
        rd = gsl_ran_beta(r, a, b);  // Radius of range size follows beta dist.
        
        // Define the grid for the species center (xc, yc)
        ix = - outer_i + 0.5 + (double)(imax + 2 * outer_i) * genrand64_real3(); // [-imax, 1+imax]
        
        // Distribution of species centers on the y-axis: 0 -> Uniform; 1 -> Normal
        if(dist_scenario == 0) {
            iy = - outer_i + 0.5 + (double)(imax + 2 * outer_i) * genrand64_real3();
        }
        else if(dist_scenario == 1) {
            variance_gauss = pow((double)imax / 5.0, 2.0);
            iy = rand_gauss((double)imax / 2.0, variance_gauss);
        }
        
        // Define the species center (x,y)
        xc = (double)ix * dx, yc = (double)iy * dx;
        
        //fprintf(fq, "%f %f\n", xc, yc);
        
        if(xc < 0.0 || xc > 1.0 || yc < 0.0 || yc > 1.0) { // Center lies outside [0,1] x [0,1]
            status = 2;
            for(i = 0; i <= imax; i++){
                p = (double)i * dx;
                dist_x0 = sqrt(pow(0.0 - xc, 2.0) + pow(p - yc, 2.0));           // Distance to the four edges
                dist_x1 = sqrt(pow(1.0 - xc, 2.0) + pow(p - yc, 2.0));
                dist_y0 = sqrt(pow(p - xc, 2.0) + pow(0.0 - yc, 2.0));
                dist_y1 = sqrt(pow(p - xc, 2.0) + pow(1.0 - yc, 2.0));
                
                if(dist_x0 < rd || dist_x1 < rd || dist_y0 < rd || dist_y1 < rd) { // Geographic range overlaps [0,1] x [0,1]
                    status = 1; break;
                }
            }
        }
        
        if(status != 2){
            count_sp++;
            count_area = 0;
            for(i = 0; i <= imax; i++){
                x = (double)i * dx;
                
                for(j = 0; j <= imax; j++){
                    y = (double)j * dx;
                    dist = sqrt(pow(x - xc, 2.0) + pow(y - yc, 2.0));
                    //printf("%f %f x_c=%f y_c=%f dist=%f\n", x, y, xc, yc, dist);
                    if(dist < rd) {
                        range_s[i][j] = 1;
                        count_area++;
                    }
                    else range_s[i][j] = 0;
                }
            }
            ptr2 = list_initial(ptr2,count_sp,imax,count_area,range_s);     // Add the species data to the list
        }
    }

    printf("tot_sp=%d, count_sp=%d; fraction of species used = %f\n", tot_sp, count_sp, (double)count_sp/(double)tot_sp );
    
    gsl_rng_free(r); // Release memory
    //fclose(fq);

}

void list_sar_ear(species *ptr2, int imax, double dx, double *u_x[], int *sar, int *ear) // Calculate SAR and EAR
{
    int i, j;
    (void)dx;
        
    for(i = 0; i <= imax; i++) {
        for(j = 0; j <= imax; j++) {
            species *q = ptr2;        // Set pointer to list head
            q = q->next;
            
            if(u_x[i][j] > EPS3){  // Cell (i,j) is converted
                //printf("u_x[%d][%d]=%f\n", i, j, u_x[i][j]);
                while(q != NULL) {
                    //printf("sp = %d sp_map[%d][%d]=%d sar=%d\n", q->species_id, i, j, q->sp_map[i][j], q->sar);
                    if(q->sp_map[i][j] == 1){ // Species is present in cell (i,j)
                        //printf("1: u_x[%d][%d]=%f sp = %d area = %d SAR:%d EAR:%d\n", i, j, u_x[i][j], q->species_id, q->area, q->sar, q->ear);
                        q->sp_map[i][j] = 0, q->area--;
                        if(q->area == 0){  // Species extinction
                            q->ear = 1;
                            ear[0] = ear[0] + 1;
                        }
                        if(q->area < (int)((double)(q->area_backup) * 0.05) && q->ear2 == 0){  // Count species with 95% range loss
                            q->ear2 = 1;
                            ear[1] = ear[1] + 1;
                        }
                        if(q->area < (int)((double)(q->area_backup) * 0.1) && q->ear3 == 0){   // Count species with 90% range loss
                            q->ear3 = 1;
                            ear[2] = ear[2] + 1;
                        }
                        if(q->sar == 0){   // First detection of the species
                            q->sar = 1;
                            *sar = *sar + 1;
                        }
                    }
                    q = q->next;
                }
            }
        }
    }
}

static int pa_point_pattern(int imax, double dx, int *u_p[], double x[], double y[])
{
    int i, j, n_pa = 0;
    
    for(i = 0; i <= imax; i++){
        for(j = 0; j <= imax; j++){
            if(u_p[i][j] == 1){
                x[n_pa] = (double)i * dx;
                y[n_pa] = (double)j * dx;
                n_pa++;
            }
        }
    }
    
    return n_pa;
}

static double translation_weight(double x1, double y1, double x2, double y2, double window_length)
{
    double overlap_x, overlap_y, overlap_area, area_window;
    
    overlap_x = window_length - fabs(x1 - x2);
    overlap_y = window_length - fabs(y1 - y2);
    if(overlap_x <= 0.0 || overlap_y <= 0.0) return 0.0;
    
    overlap_area = overlap_x * overlap_y;
    area_window = window_length * window_length;
    
    return area_window / overlap_area;
}

int pa_ripley_lr(int imax, double dx, int *u_p[], double dr, int nbin, double *lr_out)
{
    int i, j, b, n_pa;
    int max_points = (imax + 1) * (imax + 1);
    double x[max_points], y[max_points];
    double pair_count[nbin], r, dist, area_window, window_length, K, L, pi = acos(-1.0), weight;
    
    for(b = 0; b < nbin; b++) pair_count[b] = 0.0;
    
    n_pa = pa_point_pattern(imax, dx, u_p, x, y);
    
    if(n_pa < 2){
        for(b = 0; b < nbin; b++) lr_out[b] = NAN;
        return 0;
    }
    
    window_length = (double)imax * dx;
    for(i = 0; i < n_pa - 1; i++){
        for(j = i + 1; j < n_pa; j++){
            dist = sqrt((x[i] - x[j]) * (x[i] - x[j]) + (y[i] - y[j]) * (y[i] - y[j]));
            weight = translation_weight(x[i], y[i], x[j], y[j], window_length);
            if(weight == 0.0) continue;
            for(b = 0; b < nbin; b++){
                r = dr * (double)(b + 1);
                if(dist <= r) pair_count[b] += weight;
            }
        }
    }
    
    area_window = window_length * window_length;
    for(b = 0; b < nbin; b++){
        r = dr * (double)(b + 1);
        K = area_window * 2.0 * pair_count[b] / ((double)n_pa * (double)(n_pa - 1));
        L = sqrt(K / pi);
        lr_out[b] = L - r;
    }
    
    return 1;
}

int pa_paircorr_g(int imax, double dx, int *u_p[], double dr, int nbin, double *g_out)
{
    int i, j, b, n_pa;
    int max_points = (imax + 1) * (imax + 1);
    double x[max_points], y[max_points];
    double pair_count[nbin], dist, area_window, window_length, r_upper, r_lower, r_mid, pi = acos(-1.0), weight;
    
    for(b = 0; b < nbin; b++) pair_count[b] = 0.0;
    
    n_pa = pa_point_pattern(imax, dx, u_p, x, y);
    
    if(n_pa < 2){
        for(b = 0; b < nbin; b++) g_out[b] = NAN;
        return 0;
    }
    
    window_length = (double)imax * dx;
    for(i = 0; i < n_pa - 1; i++){
        for(j = i + 1; j < n_pa; j++){
            dist = sqrt((x[i] - x[j]) * (x[i] - x[j]) + (y[i] - y[j]) * (y[i] - y[j]));
            weight = translation_weight(x[i], y[i], x[j], y[j], window_length);
            if(weight == 0.0) continue;
            for(b = 0; b < nbin; b++){
                r_upper = dr * (double)(b + 1);
                r_lower = r_upper - dr;
                if(dist > r_lower && dist <= r_upper) pair_count[b] += weight;
            }
        }
    }
    
    area_window = window_length * window_length;
    for(b = 0; b < nbin; b++){
        r_upper = dr * (double)(b + 1);
        r_mid = r_upper - dr / 2.0;
        g_out[b] = area_window * pair_count[b] / (pi * r_mid * dr * (double)n_pa * (double)(n_pa - 1));
    }
    
    return 1;
}

void establish_pa_species_richness(species *ptr2, int imax, double *u_x[], int *u_p[], double max_pa_frac)
{
    int i, j, n, target_pa, current_pa, best_i, best_j, best_count;
    int species_count[imax+1][imax+1];
    species *q;
    
    target_pa = (int)ceil(max_pa_frac * pow((double)(imax + 1), 2.0));
    current_pa = 0;
    for(i = 0; i <= imax; i++) {
        for(j = 0; j <= imax; j++) {
            species_count[i][j] = 0;
            if(u_p[i][j] == 1) current_pa++;
        }
    }
    
    q = ptr2;
    q = q->next;
    while(q != NULL) {
        for(i = 0; i <= imax; i++) {
            for(j = 0; j <= imax; j++) {
                if(u_x[i][j] <= EPS3 && u_p[i][j] == 0 && q->sp_map[i][j] == 1) species_count[i][j]++;
            }
        }
        q = q->next;
    }
    
    for(n = current_pa; n < target_pa; n++) {
        best_i = -1, best_j = -1, best_count = -1;
        for(i = 0; i <= imax; i++) {
            for(j = 0; j <= imax; j++) {
                if(u_x[i][j] <= EPS3 && u_p[i][j] == 0 && species_count[i][j] > best_count) {
                    best_i = i, best_j = j, best_count = species_count[i][j];
                }
            }
        }
        if(best_i < 0 || best_j < 0) break;
        u_p[best_i][best_j] = 1;
    }
}

int judge_convergence(int count, double area_history[], double pop_history[])
{
    int i;
    double area_min, area_max, pop_min, pop_max;
    
    if(count < CONV_WINDOW_COUNT) return 0;
    
    area_min = area_max = area_history[0];
    pop_min = pop_max = pop_history[0];
    for(i = 1; i < count; i++) {
        if(area_history[i] < area_min) area_min = area_history[i];
        if(area_history[i] > area_max) area_max = area_history[i];
        if(pop_history[i] < pop_min) pop_min = pop_history[i];
        if(pop_history[i] > pop_max) pop_max = pop_history[i];
    }
    
    if(area_max - area_min < CONV_TOL_AREA && pop_max - pop_min < CONV_TOL_POP) return 1;
    else return 0;
}

void free_list(species *ptr2)
{
    species *ptr3;
    
    ptr2 = ptr2->next;
    while (ptr2 != NULL) {
        ptr3 = ptr2->next;
        free(ptr2);
        ptr2 = ptr3;
    }
}

void data_connection(int **data[], double *data_time[], int *cover[], int scenario, int a_step, double dw2, int TRIAL, double max_pa_frac, double lu_level_pa_est){
    
    FILE *fq4, *fq5, *fq6;
    fq4 = fopen("data_result.dat","w"), fq5 = fopen("data_timing.dat","w"), fq6 = fopen("data_coverage.dat","w");
    
    int i, j, k;
    double increment, increment_2, increment_1, increment0, increment1, increment2, ma, coverage, max_x_cover;
    
    for(i = 0; i < scenario; i++){
        for(k = 0; k < 4; k++){
            for(j = 1; j < a_step; j++){
                if(data[i][j][k] < data[i][j-1][k]){
                    data[i][j][k] = data[i][j-1][k];
                }
            }
        }
    }
    
    fprintf(fq4, "### Max PA fraction=%.5f PA is established after LU fraction = %.5f\n", max_pa_frac, lu_level_pa_est);
    for(i = 0; i < scenario; i++){
        for(k = 0; k < 4; k++){ // 100% extinction EAR (k=0), 95% EAR (k=1), 90% EAR (k=2), and SAR (k=3)
            fprintf(fq4, "### SCENARIO=%d k=%d\n0.000000 0 0.000000 0.000000 nan\n", i, k);
            for(j = 0; j < a_step - 1; j++){
                if(j > 2 && j < a_step - 3) {
                    increment_2 = ((double)data[i][j-2][k] - (double)data[i][j-3][k]) / (double)a_step;
                    increment_1 = ((double)data[i][j-1][k] - (double)data[i][j-2][k]) / (double)a_step;
                    increment0 = ((double)data[i][j][k] - (double)data[i][j-1][k]) / (double)a_step;
                    increment1 = ((double)data[i][j+1][k] - (double)data[i][j][k]) / (double)a_step;
                    increment2 = ((double)data[i][j+2][k] - (double)data[i][j+1][k]) / (double)a_step;
                    ma = (increment_2 + increment_1 + increment0 + increment1 + increment2) / 5.0; // Five-point moving average
                }
                else ma = NAN;
                if(j > 0) increment = ((double)data[i][j][k] - (double)data[i][j-1][k]) / (double)a_step;
                else increment = 0.0;
                fprintf(fq4, "%f %d %f %f %f\n", dw2 * (double)(j+1), data[i][j][k], (double)data[i][j][k] / (double)TRIAL, increment, ma);
            }
            fprintf(fq4, "\n\n");
            if(k == 0){
                for(j = 1; j < 20; j++){
                    fprintf(fq5, "%f %f\n", 0.05 * (double)j, (double)data_time[i][j] / (double)TRIAL);
                }
                fprintf(fq5, "\n\n");
            }
        }
        max_x_cover = 0.0;
        fprintf(fq6, "### SCENARIO=%d COVERAGE_LEVEL=%.3f\n", i, COVERAGE_LEVEL);
        for(j = 0; j < a_step - 1; j++){
            coverage = (double)cover[i][j] / (double)TRIAL;
            if(coverage >= COVERAGE_LEVEL) max_x_cover = dw2 * (double)(j+1);
            fprintf(fq6, "%f %f %d\n", dw2 * (double)(j+1), coverage, cover[i][j]);
        }
        fprintf(fq6, "#MAX_X_COVER=%f\n\n", max_x_cover);
    }
    
    fclose(fq4), fclose(fq5), fclose(fq6);
}

void ripley_connection(double *ripley_sum[], int ripley_count[], int scenario, int nbin, double dr)
{
    FILE *fq;
    int i, j;
    double mean_lr;
    
    fq = fopen("data_ripley.dat","w");
    
    for(i = 0; i < scenario; i++){
        fprintf(fq, "### SCENARIO=%d VALID_RUNS=%d\n", i, ripley_count[i]);
        for(j = 0; j < nbin; j++){
            if(ripley_count[i] > 0) mean_lr = ripley_sum[i][j] / (double)ripley_count[i];
            else mean_lr = NAN;
            fprintf(fq, "%f %f\n", dr * (double)(j + 1), mean_lr);
        }
        fprintf(fq, "\n\n");
    }
    
    fclose(fq);
}

void paircorr_connection(double *paircorr_sum[], int paircorr_count[], int scenario, int nbin, double dr)
{
    FILE *fq;
    int i, j;
    double mean_g;
    
    fq = fopen("data_paircorr.dat","w");
    
    for(i = 0; i < scenario; i++){
        fprintf(fq, "### SCENARIO=%d VALID_RUNS=%d\n", i, paircorr_count[i]);
        for(j = 0; j < nbin; j++){
            if(paircorr_count[i] > 0) mean_g = paircorr_sum[i][j] / (double)paircorr_count[i];
            else mean_g = NAN;
            fprintf(fq, "%f %f\n", dr * (double)(j + 1), mean_g);
        }
        fprintf(fq, "\n\n");
    }
    
    fclose(fq);
}

void PDE(species *ptr2, int tmax, int imax, double dt, double Dx, double Dy, double r_x, double K, double dx, double dx2, double *u_x[], int *u_p[], double **u_matrix[], int PRINT, int step_lu, int step_pa, int **data[], double *data_time[], int *cover[], int scenario, int a_step, double lu_level_pa_est, double max_pa_frac){
    
    FILE *fp, *fp2;
    fp = fopen("data_urban_map.dat","w"), fp2 = fopen("data_ear.dat","a");
    
    fprintf(fp2, "#STEP urbanization=%d\n", step_lu);
    
    int i, i2, j, t, order, ix, iy, pa_established, conv_sample_step, conv_count, conv_index;
    double area_lu = 0.0, area_pa = 0.0, population_x, final_cover_area;
    double area_history[CONV_WINDOW_COUNT], pop_history[CONV_WINDOW_COUNT];
    
    int data_temp[a_step][4]; // Keep this on the stack; do not replace with calloc here
    for(i = 0; i < a_step; i++ ){
        for(j = 0; j < 4; j++){
            data_temp[i][j] = 0;
        }
    }
    
    int sar = 0;
    int ear[3];
    
    i2 = 1, ear[0] = ear[1] = ear[2] = 0;
    pa_established = 0;
    conv_sample_step = (int)(CONV_SAMPLE_TIME / dt + 0.5);
    if(conv_sample_step < 1) conv_sample_step = 1;
    conv_count = 0, conv_index = 0;
    
    // Define the initial land-use center (xc, yc)
    ix = 0.5 + imax * genrand64_real3(); // [0, imax]
    iy = 0.5 + imax * genrand64_real3();
    u_x[ix][iy] = 10.0;
    urban_pa_area(imax, u_x, u_p, &area_lu, &area_pa);
    
    for(t = 1; t < tmax; t++){ // Time loop
        if(t % step_lu == 0 && (area_lu + area_pa) < 0.99){ // A new land-use event occurs
            do{
                // Define the land-use center (xc, yc)
                ix = 0.5 + imax * genrand64_real3(); // [0, imax]
                iy = 0.5 + imax * genrand64_real3();
            } while(u_x[ix][iy] > EPS3 || u_p[ix][iy] == 1);             // Place the new land-use seed in a non-converted, non-protected cell
            u_x[ix][iy] = K;
        }
        // Introduce new protected cells once land use exceeds "lu_level_pa_est", up to "max_pa_frac"
        if(t % step_pa == 0 && pa_established == 0 && area_lu > lu_level_pa_est && area_pa < max_pa_frac && (area_lu + area_pa) < 0.99){
            establish_pa_species_richness(ptr2, imax, u_x, u_p, max_pa_frac);
            pa_established = 1;
        }
        
        // RK loop for fourth-order accuracy
        for(order = 1; order <= 4; order++){
            
            #pragma omp parallel for num_threads(12) private (j)
            for(i = 1; i < imax; i++){
                for(j = 1; j < imax; j++){
                    if(u_p[i][j] == 0){ // Location (i,j) is not protected
                        // Runge-Kutta method
                        u_matrix[i][j][order] = dt * rk_2d(i,j,order,imax,Dx,Dy,r_x,K,dx,u_x,u_matrix);
                        //printf( "RUNGE=KUTTAt=%d %d %d %d %f \n", t, i, j, order, u_matrix[i][j][order] );
                    }
                }
            }
        }
        
        // Update u_x
        for(i = 1; i < imax; i++) {
            for(j = 1; j < imax; j++) {
                // Runge-Kutta method
                if(u_p[i][j] == 0){
                    u_x[i][j] = u_x[i][j] + (u_matrix[i][j][1] + 2.0 * ( u_matrix[i][j][2] + u_matrix[i][j][3] ) + u_matrix[i][j][4]) / 6.0;
                }
            }
        }
        
        // Boundary conditions
        for (i = 1; i < imax; i++){
            // Neumann
            u_x[imax][i] = u_x[imax-2][i];
            u_x[i][imax] = u_x[i][imax-2];
            u_x[i][0] = u_x[i][2];
            u_x[0][i] = u_x[2][i];
        }
        
        // Four corners
        u_x[0][0] = (u_x[2][0] + u_x[0][2]) / 2.0;
        u_x[0][imax] = (u_x[0][imax-2] + u_x[2][imax]) / 2.0;
        u_x[imax][0] =  (u_x[imax-2][0] + u_x[imax][2]) / 2.0;
        u_x[imax][imax] = (u_x[imax-2][imax] + u_x[imax][imax-2]) / 2.0;
        
        population_x = abundance_2d(u_x,imax,dx);
        //printf( "t=%d x = %f dif=%f\n", t, population_x, fabs(population_x - population_x_prev) );
        
        // Check convergence
        if(population_x < EPS2) {
            printf("Extinction!\n");
            break;
        }
                
        list_sar_ear(ptr2,imax,dx,u_x,&sar,ear);
        urban_pa_area(imax,u_x,u_p,&area_lu,&area_pa); // Calculate the fractions of converted and protected area
             
        if(t % 100 == 0) printf("t=%f area_lu=%f area_pa=%f population_x=%f\n", (double)t * dt, area_lu, area_pa, population_x);
        
        // Pool EAR/SAR data over the trial loop by transformed-area bin; average at the end
        for(i = 0; i < a_step; i++){
            if(dx2 * (double)i <= area_lu && dx2 * (double)(i + 1) > area_lu){
                data_temp[i][0] = ear[0];
                data_temp[i][1] = ear[1];
                data_temp[i][2] = ear[2];
                data_temp[i][3] = sar;      // SAR
            }
        }
        
        //if(area_lu + area_pa >= 0.99) break;
          
        if(area_lu + area_pa >= 0.99 && (double)t * dt > 8.0) break;
        
        if(t % conv_sample_step == 0){
            area_history[conv_index] = area_lu;
            pop_history[conv_index] = population_x;
            if(conv_count < CONV_WINDOW_COUNT) conv_count++;
            conv_index = (conv_index + 1) % CONV_WINDOW_COUNT;
            
            if((double)t * dt > CONV_MIN_TIME && judge_convergence(conv_count, area_history, pop_history) == 1) break;
        }
          
        // Input timing data
        if(area_lu >= 0.04999 * (double)i2) {
            data_time[scenario][i2] += (double)t * dt;
            i2++;
        }
        
        //fprintf(fp2, "%f %f %f %d %d\n", (double)t * dt, population_x, area, sar, ear);
    
    } // End of the time loop
    
    //fprintf(fp2, "\n\n");
    //printf("SCENARIO=%d t=%d Area = %f SAR=%d EAR=%d\n", scenario, t, area, sar, ear);
    
    for(i = 0; i < a_step; i++){ // Write output data
        data[scenario][i][0] += data_temp[i][0];
        data[scenario][i][1] += data_temp[i][1];
        data[scenario][i][2] += data_temp[i][2];
        data[scenario][i][3] += data_temp[i][3];
    }
    
    final_cover_area = area_lu + area_pa;
    for(i = 0; i < a_step - 1; i++){
        if(final_cover_area >= dx2 * (double)(i + 1)) cover[scenario][i]++;
    }
    
    // Print data
    if(PRINT == 1){
        for(i = 0; i <= imax; i++) {
            for(j = 0; j <= imax; j++) {
                fprintf(fp, "%f %f %f %d\n", (double)i * dx, (double)j * dx, u_x[i][j], u_p[i][j]);
            }
            fprintf(fp, "\n");
        }fprintf(fp, "\n");
    }
    
    fclose(fp), fclose(fp2);
    
}

double rk_2d(int l, int k, int rk, int lmax, double Dx, double Dy, double r, double K, double dx, double *x[], double **matrix_x[]){
    
    int i, j;
    double value, x2[lmax+1][lmax+1];
        
    if (rk == 1){                      // First RK stage
        for (i = 0; i <= lmax; i++){
            for (j = 0; j <= lmax; j++){
                x2[i][j] = x[i][j];
            }
        }
    }
    else if (rk == 2 || rk == 3){
        for (i = 0; i <= lmax; i++){
            for (j = 0; j <= lmax; j++){
                x2[i][j] = x[i][j] + matrix_x[i][j][rk-1] / 2.0;
            }
        }
    }
    else if (rk == 4){
        for (i = 0; i <= lmax; i++){
            for (j = 0; j <= lmax; j++){
                x2[i][j] = x[i][j] + matrix_x[i][j][3];
            }
        }
    }
    
    value = Dx * (x2[l+1][k] - 2.0 * x2[l][k] + x2[l-1][k]) / (dx * dx) + Dy * (x2[l][k+1] - 2.0 * x2[l][k] + x2[l][k-1]) / (dx * dx) + r * x2[l][k] * (1.0 - x2[l][k] / K);

    return value;
    
}

void urban_pa_area(int imax, double *u_x[], int *u_p[], double *area_lu, double *area_pa){
    
    int i, j;
    int count_urban_area = 0;
    int count_pa = 0;
    
    for(i = 0; i <= imax; i++) {
        for(j = 0; j <= imax; j++) {
            if(u_x[i][j] > EPS3) count_urban_area++;
            if(u_p[i][j] == 1) count_pa++;
        }
    }

    *area_lu = (double)count_urban_area / pow((double)(imax + 1), 2.0);
    *area_pa = (double)count_pa / pow((double)(imax + 1), 2.0);
}

double abundance_2d(double *u[], int xmax, double dx){
    
    int i, j;
    double u_sum = 0.0;
    
    for( i = 0; i <= xmax; i++ ){
        for( j = 0; j <= xmax; j++ ){
            u_sum += u[i][j] * dx * dx;
        }
    }
    
    return u_sum;
}


void release_memory_int_2d(int *u[], int imax){
    
    int i;
    
    // Memory release
    for(i = 0; i < imax; i++) free(u[i]);
    
    free(u);
}

void release_memory_int_3d(int **u[], int imax, int jmax){
    
    int i, j;
    
    // Memory release
    for(i = 0; i < imax; i++){
        for(j = 0; j < jmax; j++){
            free(u[i][j]);
        }
    }
    
    free(u);
    
}

void release_memory_double_2d( double *u[], int imax ){
    
    int i;
    
    // Memory release
    for( i = 0; i < imax; i++ ){
        free(u[i]);
    }
    
    free(u);
}

void release_memory_double_3d( double **u[], int imax, int jmax){
    
    int i, j;
    
    // Memory release
    for( i = 0; i < imax; i++ ){
        for( j = 0; j < jmax; j++ ){
            free(u[i][j]);
        }
    }
    
    free(u);
    
}

int rand_po(double u) // Poisson generator; Devroye, Luc (1986). "Discrete Univariate Distributions" p. 504
{
    
    int k;
    double p, s, rand;
    
    k = 0;
    p = exp( -u );
    s = p;
    rand = genrand64_real1();

    while( rand > s ){
        k++;
        p = p * u / (double)k;
        s = s + p;
    }
    
    return k;
    
}

double rand_gauss(double mu, double variance){
    
    double x1, x2, m_n1;
    
    x1 = genrand64_real3();
    x2 = genrand64_real3();
    
    // Box-Muller's method; generate p(x)~N(0, 1), see https://en.wikipedia.org/wiki/Box%E2%80%93Muller_transform
    m_n1 = sqrt( - 2.0 * log(x1) ) * sin( 2.0 * M_PI * x2 );
    return m_n1 * sqrt( variance ) + mu; // Draw u ~ N(mu, variance)
    
}

/*
 A 64-bit Mersenne Twister generator: http://www.math.hiroshima-u.ac.jp/~m-mat/MT/emt.html
 */

#define NN 312
#define MM 156
#define MATRIX_A 0xB5026F5AA96619E9ULL
#define UM 0xFFFFFFFF80000000ULL /* Most significant 33 bits */
#define LM 0x7FFFFFFFULL /* Least significant 31 bits */

/* The array for the state vector */
static unsigned long long mt[NN];
/* mti==NN+1 means mt[NN] is not initialized */
static int mti=NN+1;

/* initializes mt[NN] with a seed */
void init_genrand64(unsigned long long seed)
{
    mt[0] = seed;
    for (mti=1; mti<NN; mti++)
        mt[mti] =  (6364136223846793005ULL * (mt[mti-1] ^ (mt[mti-1] >> 62)) + mti);
}

/* generates a random number on [0, 2^64-1]-interval */
unsigned long long genrand64_int64(void)
{
    int i;
    unsigned long long x;
    static unsigned long long mag01[2]={0ULL, MATRIX_A};
    
    if (mti >= NN) { /* generate NN words at one time */
        
        /* if init_genrand64() has not been called, */
        /* a default initial seed is used     */
        if (mti == NN+1)
            init_genrand64(5489ULL);
        
        for (i=0;i<NN-MM;i++) {
            x = (mt[i]&UM)|(mt[i+1]&LM);
            mt[i] = mt[i+MM] ^ (x>>1) ^ mag01[(int)(x&1ULL)];
        }
        for (;i<NN-1;i++) {
            x = (mt[i]&UM)|(mt[i+1]&LM);
            mt[i] = mt[i+(MM-NN)] ^ (x>>1) ^ mag01[(int)(x&1ULL)];
        }
        x = (mt[NN-1]&UM)|(mt[0]&LM);
        mt[NN-1] = mt[MM-1] ^ (x>>1) ^ mag01[(int)(x&1ULL)];
        
        mti = 0;
    }
    
    x = mt[mti++];
    
    x ^= (x >> 29) & 0x5555555555555555ULL;
    x ^= (x << 17) & 0x71D67FFFEDA60000ULL;
    x ^= (x << 37) & 0xFFF7EEE000000000ULL;
    x ^= (x >> 43);
    
    return x;
}

/* generates a random number on [0, 2^63-1]-interval */
long long genrand64_int63(void)
{
    return (long long)(genrand64_int64() >> 1);
}

/* generates a random number on [0,1]-real-interval */
double genrand64_real1(void)
{
    return (genrand64_int64() >> 11) * (1.0/9007199254740991.0);
}

/* generates a random number on [0,1)-real-interval */
double genrand64_real2(void)
{
    return (genrand64_int64() >> 11) * (1.0/9007199254740992.0);
}

/* generates a random number on (0,1)-real-interval */
double genrand64_real3(void)
{
    return ((genrand64_int64() >> 12) + 0.5) * (1.0/4503599627370496.0);
}
