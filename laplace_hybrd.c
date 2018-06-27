/****************************************************************
 * Laplace MPI Template C Version
 *
 * T is initially 0.0
 * Boundaries are as follows
 *
 *                T                      4 sub-grids
 *   0  +-------------------+  0    +-------------------+
 *      |                   |       |                   |
 *      |                   |       |-------------------|
 *      |                   |       |                   |
 *   T  |                   |  T    |-------------------|
 *      |                   |       |                   |
 *      |                   |       |-------------------|
 *      |                   |       |                   |
 *   0  +-------------------+ 100   +-------------------+
 *      0         T       100
 *
 * Each PE only has a local subgrid.
 * Each PE works on a sub grid and then sends
 * Its boundaries to neighbors.
 *
 *  John Urbanic, PSC 2014
 *
 *******************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <sys/time.h>
#include <mpi.h>

#define COLUMNS_GLOBAL   10000                  // this is a "global" column count
#define ROWS_GLOBAL      10000                  // this is a "global" row count
//#define NPES                4                  // number of processors
//#define ROWS             (ROWS_GLOBAL/NPES)    // number of real local rows
//#define COLUMNS          (COLUMNS_GLOBAL/NPES) // number of real local columns

// communication tags
#define DOWN             100
#define UP               101
#define RIGHT            102
#define LEFT             103

#define MAX_TEMP_ERROR   0.01

void initialize(int npes, int my_PE_num, int PEi, int PEj, int ROWS, int COLUMNS, double** Temperature);
void track_progress(int iter, int ROWS, int COLUMNS, double** Temperature);

int main(int argc, char *argv[]) {

    int i, j;
    int max_iterations;
    int iteration=1;
    double dt;
    struct timeval start_time, stop_time, elapsed_time;

    int        nnpes; // = NPES*NPES;
    int        npes; // = NPES*NPES;
    int        my_PE_num;           // my PE number
    int        PEi, PEj;
    int        ROWS;
    int        COLUMNS;
    double     **Temperature;
    double     **Temperature_last;
    double     *temp_to_send; 
    double     *temp_to_recv;
    int        PE_UP, PE_DOWN, PE_RIGHT, PE_LEFT;

    double     dt_global = 100;       // delta t across all PEs
    MPI_Status status;              // status returned by MPI calls

    // the usual MPI startup routines
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_PE_num);
    MPI_Comm_size(MPI_COMM_WORLD, &nnpes);

    // decompose the proc num and compute the size of a subdomain.
    // Assumed configuration:
    //    0             1               ... npes-1
    //    npes          npes+1          ... 2*npes-1
    //    ...           ...             ... ...
    //    (npes-1)*npes (npes-1)*npes+1 ... npes*npes-1
    npes = sqrt( nnpes );
    PEj = my_PE_num % npes;
    PEi = (my_PE_num - PEj) / npes;
    ROWS = ROWS_GLOBAL / npes;
    COLUMNS = COLUMNS_GLOBAL / npes;

    // compute neighbor proc numbers.
    PE_LEFT  = (PEj-1) + (npes*PEi);
    PE_RIGHT = (PEj+1) + (npes*PEi);
    PE_DOWN  = PEj + npes*(PEi+1);
    PE_UP    = PEj + npes*(PEi-1);   

    /* allocate the array */
    Temperature = malloc((ROWS+2) * sizeof(*Temperature));
    Temperature_last = malloc((ROWS+2) * sizeof(*Temperature_last));
    for (i=0; i<ROWS+2; i++)
    {
        Temperature[i] = malloc((COLUMNS+2) * sizeof(*Temperature[i]));
        Temperature_last[i] = malloc((COLUMNS+2) * sizeof(*Temperature_last[i]));
    }
    temp_to_send = malloc( ROWS * sizeof(double) );
    temp_to_recv = malloc( ROWS * sizeof(double) );

    // verify only NPES PEs are being used
    /*
    if(npes != NPES) {
      if(my_PE_num==0) {
        printf("This code must be run with %d PEs\n", NPES);
      }
      MPI_Finalize();
      exit(1);
    }
    */

    // PE 0 asks for input

    if (my_PE_num == 0) {
      printf("How many iterations?\n");
      //fflush(stdout); // Not always necessary, but can be helpful // Idk what is it for
      scanf( "%d", &max_iterations);
    }

    // bcast max iterations to other PEs
    MPI_Bcast(&max_iterations, 1, MPI_INT, 0, MPI_COMM_WORLD);


    if (my_PE_num==0) gettimeofday(&start_time,NULL);

    // Initialize Temperature array with boundary conditions for each PE
    initialize(npes, my_PE_num, PEi, PEj, ROWS, COLUMNS, Temperature_last);

   #pragma acc data copy(Temperature_last), create(Temperature)

    while ( dt_global > MAX_TEMP_ERROR && iteration <= max_iterations ) {

        #pragma acc kernels

        // main calculation: average my four neighbors
        for(i = 1; i <= ROWS; i++) {
            for(j = 1; j <= COLUMNS; j++) {
                Temperature[i][j] = 0.25 * (Temperature_last[i+1][j] + Temperature_last[i-1][j] +
                                            Temperature_last[i][j+1] + Temperature_last[i][j-1]);
            }
        }

        #pragma acc update host(Temperature[1:1][1:COLUMNS], Temperature[ROWS:1][1:COLUMNS]) 

        if ( PEi != npes-1 ){
            MPI_Send(&Temperature[ROWS][1], COLUMNS, MPI_DOUBLE, PE_DOWN, DOWN, MPI_COMM_WORLD);
            //printf("Sending my_PE_num = %d PE_DOWN = %d\n", my_PE_num, PE_DOWN);
        }

        if ( PEi != 0 ){
            MPI_Recv(&Temperature_last[0][1], COLUMNS, MPI_DOUBLE, PE_UP, DOWN, MPI_COMM_WORLD, &status);
            //printf("Receiving my_PE_num = %d PE_UP = %d\n", my_PE_num, PE_UP);
        }

        if ( PEi != 0 ){
            MPI_Send(&Temperature[1][1], COLUMNS, MPI_DOUBLE, PE_UP, UP, MPI_COMM_WORLD);
            //printf("Sending my_PE_num = %d PE_UP = %d\n", my_PE_num, PE_UP);
        }

        if ( PEi != npes-1 ){
            MPI_Recv(&Temperature_last[ROWS+1][1], COLUMNS, MPI_DOUBLE, PE_DOWN, UP, MPI_COMM_WORLD, &status);
            //printf("Receiving my_PE_num = %d PE_DOWN = %d\n", my_PE_num, PE_DOWN);
        }
     
        #pragma acc update device(Temperature_last[0:1][1:COLUMNS], Temperature_last[ROWS+1:1][1:COLUMNS]) 

        //MPI_Barrier(MPI_COMM_WORLD);

        #pragma acc update host(Temperature[ROWS:1][ROWS+1:1], Temperature[ROWS:1][ROWS+1:COLUMNS])

        if ( PEj != npes-1 ){
            for(i = 0; i < ROWS; i++ )
               temp_to_send[i] = Temperature[i+1][COLUMNS];
            MPI_Send(temp_to_send, ROWS, MPI_DOUBLE, PE_RIGHT, RIGHT, MPI_COMM_WORLD);
            //printf("Sending my_PE_num = %d PE_RIGHT = %d\n", my_PE_num, PE_RIGHT);
        }

        if ( PEj != 0 ){
            MPI_Recv(temp_to_recv, ROWS, MPI_DOUBLE, PE_LEFT, RIGHT, MPI_COMM_WORLD, &status);
            for(i = 0; i < ROWS; i++ )
               Temperature_last[i+1][0] = temp_to_recv[i];
            //printf("Receiving my_PE_num = %d PE_LEFT = %d\n", my_PE_num, PE_LEFT);
        }

        if ( PEj != 0 ){
            for(i = 0; i < ROWS; i++ )
               temp_to_send[i] = Temperature[i+1][1];
            MPI_Send(temp_to_send, ROWS, MPI_DOUBLE, PE_LEFT, LEFT, MPI_COMM_WORLD);
            //printf("Sending my_PE_num = %d PE_LEFT = %d\n", my_PE_num, PE_LEFT);
        }

        if ( PEj != npes-1 ){
            MPI_Recv(temp_to_recv, ROWS, MPI_DOUBLE, PE_RIGHT, LEFT, MPI_COMM_WORLD, &status);
            for(i = 0; i < ROWS; i++ )
               Temperature_last[i+1][COLUMNS+1] = temp_to_recv[i];
            //printf("Receiving my_PE_num = %d PE_RIGHT = %d\n", my_PE_num, PE_RIGHT);
        }

        #pragma acc update device(Temperature_last[ROWS:1][ROWS+1:0], Temperature_last[ROWS:1][ROWS+1:COLUMNS+1])

        //printf("my_PE_num = %d\n", my_PE_num);
        //MPI_Barrier(MPI_COMM_WORLD);

        dt = 0.0;

        #pragma acc kernels

        for(i = 1; i <= ROWS; i++){
            for(j = 1; j <= COLUMNS; j++){
                dt = fmax( fabs(Temperature[i][j]-Temperature_last[i][j]), dt);
                Temperature_last[i][j] = Temperature[i][j];
            }
        }

        // find global dt
        MPI_Reduce(&dt, &dt_global, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
        MPI_Bcast(&dt_global, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        // periodically print test values - only for PE in lower corner
        if((iteration % 100) == 0) {
            if (my_PE_num == npes-1){
                track_progress(iteration,ROWS, COLUMNS, Temperature_last);
            }
        }

        iteration++;
    }

    // Slightly more accurate timing and cleaner output
    MPI_Barrier(MPI_COMM_WORLD);

    // PE 0 finish timing and output values
    if (my_PE_num==0){
        gettimeofday(&stop_time,NULL);
        timersub(&stop_time, &start_time, &elapsed_time);

        printf("\nMax error at iteration %d was %f\n", iteration-1, dt_global);
        printf("Total time was %f seconds.\n", elapsed_time.tv_sec+elapsed_time.tv_usec/1000000.0);
    }

    // print temperature at the point (7500,9950)
     if (my_PE_num==79){
        printf("Global coord [7500,9950] is %f \n:", Temperature[500][950]);
     }

    // Free up memory allocated to temperature arrays.
    for (i=0; i<ROWS+2; i++)
    {
        free(Temperature[i]);
        free(Temperature_last[i]);
    }
    free(Temperature);
    free(Temperature_last);
    free(temp_to_send);
    free(temp_to_recv);

    MPI_Finalize();
}



void initialize(int npes, int my_PE_num, int PEi, int PEj, int ROWS, int COLUMNS, double** Temperature_last){

    double tMin_b, tMax_b, tMin_l, tMax_l;  //Local boundary limits
    int i,j;
    
    for(i = 0; i <= ROWS+1; i++){
        for (j = 0; j <= COLUMNS+1; j++){
            Temperature_last[i][j] = 0.0;
        }
    }

    // Local boundry condition endpoints
    /*
    tMin = (my_PE_num)*100.0/npes;
    tMax = (my_PE_num+1)*100.0/npes;
    */

    /*    
    double re_my_PE_num;
    if (my_PE_num < nnpes)
        re_my_PE_num = 0
    if (nnpes <= my_PE_num < 2*nnpes)
        re_my_PE_num = 1
    if (2*nnpes <= my_PE_num < 3*nnpes)
        re_my_PE_num = 2
    */
     
    tMin_b = (PEj * 100.0) / npes;
    tMax_b = ((PEj + 1) * 100.0) / npes;

    tMin_l = (PEi * 100.0) / npes;
    tMax_l = ((PEi + 1) * 100.0) / npes;

    // Left and right boundaries
    for (i = 0; i <= ROWS+1; i++) {
      Temperature_last[i][0] = 0.0;
      Temperature_last[i][COLUMNS+1] = tMin_l + ((tMax_l-tMin_l)/(ROWS+1))*i;
    }

    // Top boundary (PE 0 only)
    if (PEi == 0)
      for (j = 0; j <= COLUMNS+1; j++)
        Temperature_last[0][j] = 0.0;

    // Bottom boundary (Last PE only)
    if (PEi == npes-1)
      for (j=0; j<=COLUMNS+1; j++)
        Temperature_last[ROWS+1][j] = tMin_b + ((tMax_b-tMin_b)/(COLUMNS+1))*j;
                                      //(((100.0/(COLUMNS+1)) * j)/npes)*(PEj+1);

//    for( i = 0; i <= ROWS+1; i++ )
//      for( j = 0; j <= COLUMNS+1; j++ )
//        fprintf(stderr, "me: %d T=%e (i,j)=%d %d\n", my_PE_num, Temperature_last[i][j], i, j);



//     if (PEi ==npes -1)
//         for (j = 0; j <= COLUMNS+1; j++)
//            fprintf(stderr, "me: %d T = %e (i,j) = %d %d\n", my_PE_num, Temperature_last[ROWS+1][j], ROWS+1, j);
     

//     if (PEj ==npes -1)
//         for (i= 0; i <= ROWS+1; i++)
//           fprintf(stderr, "me: %d T = %e (i,j) = %d %d\n", my_PE_num, Temperature_last[i][COLUMNS+1], i, COLUMNS+1);


//     fprintf(stderr,"my_PE_num: %d PEi = %d PEj = %d\n", my_PE_num, PEi, PEj); 

}

// only called by last PE
void track_progress(int iteration, int ROWS, int COLUMNS, double** Temperature) {

    int i;

    printf("---------- Iteration number: %d ------------\n", iteration);

    // output global coordinates so user doesn't have to understand decomposition
    for(i = 5; i >= 0; i--) {
      printf("[%d,%d]: %5.2f  ", ROWS_GLOBAL-i, COLUMNS-i, Temperature[ROWS-i][COLUMNS-i]);
    }
    printf("\n");
}
