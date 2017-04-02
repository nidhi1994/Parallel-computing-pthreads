/* Gaussian elimination code.
 * Author: Naga Kandasamy
 * Date created: 02/07/2014
 * Date of last update: 01/30/2017
 * Compile as follows: gcc -o gauss_eliminate gauss_eliminate.c compute_gold.c -lpthread -std=c99 -lm
 */

#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <sys/time.h>
#include <string.h>
#include <math.h>
#include <pthread.h>
#include "gauss_eliminate.h"

#define MIN_NUMBER 2
#define MAX_NUMBER 50
#define NUM_THREADS 32

typedef struct args_for_thread_s{
      int thread_id; // The thread ID
      Matrix *U;
      int cnt;
      int k;
} ARGS_FOR_THREAD;

/* Function prototypes. */
extern int compute_gold (float *, unsigned int);
Matrix allocate_matrix (int num_rows, int num_columns, int init);
void gauss_eliminate_using_pthreads (Matrix);
int perform_simple_check (const Matrix);
void *division(void*);
void *elimination(void*);
void print_matrix (const Matrix);
float get_random_number (int, int);
int check_results (float *, float *, unsigned int, float);


int
main (int argc, char **argv)
{
  /* Check command line arguments. */
  if (argc > 1)
    {
      printf ("Error. This program accepts no arguments. \n");
      exit (0);
    }
    float ts,tp;
  /* Matrices for the program. */
  Matrix A;     // The input matrix
  Matrix U_reference;   // The upper triangular matrix computed by the reference code
  Matrix U_mt;      // The upper triangular matrix computed by the pthread code

  /* Initialize the random number generator with a seed value. */
  srand (time (NULL));

  /* Allocate memory and initialize the matrices. */
  A = allocate_matrix (MATRIX_SIZE, MATRIX_SIZE, 1);  // Allocate and populate a random square matrix
  U_reference = allocate_matrix (MATRIX_SIZE, MATRIX_SIZE, 0);  // Allocate space for the reference result
  U_mt = allocate_matrix (MATRIX_SIZE, MATRIX_SIZE, 0); // Allocate space for the multi-threaded result

//printf("original\n");
  /* Copy the contents of the A matrix into the U matrices. */
  for (int i = 0; i < A.num_rows; i++)
    {
      for (int j = 0; j < A.num_rows; j++)
  {
    U_reference.elements[A.num_rows * i + j] = A.elements[A.num_rows * i + j];
    U_mt.elements[A.num_rows * i + j] = A.elements[A.num_rows * i + j];
   // printf("  %f  ", U_mt.elements[A.num_rows * i + j]);
  }
    }
    U_mt.num_rows=A.num_rows;

  printf ("Performing gaussian elimination using the reference code. \n");
  struct timeval start, stop;
  gettimeofday (&start, NULL);
  int status = compute_gold (U_reference.elements, A.num_rows);
  gettimeofday (&stop, NULL);
  printf ("CPU run time = %0.2f s. \n",
    (float) (stop.tv_sec - start.tv_sec +
       (stop.tv_usec - start.tv_usec) / (float) 1000000));
ts=(float) (stop.tv_sec - start.tv_sec +
       (stop.tv_usec - start.tv_usec) / (float) 1000000);
  if (status == 0)
    {
      printf("Failed to convert given matrix to upper triangular. Try again. Exiting. \n");
      exit (0);
    }
  status = perform_simple_check (U_reference);  // Check that the principal diagonal elements are 1 
  if (status == 0)
    {
      printf ("The upper triangular matrix is incorrect. Exiting. \n");
      exit (0);
    }
  printf ("Single-threaded Gaussian elimination was successful. \n");

  /* Perform the Gaussian elimination using pthreads. The resulting upper triangular matrix should be returned in U_mt */
  gettimeofday (&start, NULL);
  gauss_eliminate_using_pthreads (U_mt);
  gettimeofday (&stop, NULL);
  printf ("CPU run time = %0.2f s. \n",
    (float) (stop.tv_sec - start.tv_sec +
       (stop.tv_usec - start.tv_usec) / (float) 1000000));
tp= (float) (stop.tv_sec - start.tv_sec +
       (stop.tv_usec - start.tv_usec) / (float) 1000000);
  /* check if the pthread result is equivalent to the expected solution within a specified tolerance. */
  int size = MATRIX_SIZE * MATRIX_SIZE;
  int res = check_results (U_reference.elements, U_mt.elements, size, 0.001f);
  printf("Speedup  %f\n", ts/tp );
  printf ("Test %s\n", (1 == res) ? "PASSED" : "FAILED");

  /* Free memory allocated for the matrices. */
  free (A.elements);
  free (U_reference.elements);
  free (U_mt.elements);

  return 0;
}


/* Write code to perform gaussian elimination using pthreads. */
void
gauss_eliminate_using_pthreads (Matrix U)
{
   pthread_t thread_id[NUM_THREADS]; // Data structure to store the thread IDs
  pthread_attr_t attributes; // Thread attributes
   pthread_attr_init(&attributes); // Initialize the thread attributes to the default values
   int i,m,n,k;
    double sum = 0; 
    ARGS_FOR_THREAD *args_for_thread[NUM_THREADS];
    int chunk_size = (int)floor((float)(U.num_rows)/(float)NUM_THREADS); // Compute the chunk size
    int num_elements=U.num_rows;
  for (int k = 0; k < num_elements; k++)
  {
  // int k=0;
   /* ARGS_FOR_THREAD *args_for_thread[NUM_THREADS];
     for(i = 0; i < NUM_THREADS; i++){
        args_for_thread[i] = (ARGS_FOR_THREAD *)malloc(sizeof(ARGS_FOR_THREAD));
        args_for_thread[i]->thread_id = i; // Provide thread ID  
        args_for_thread[i]->U=&U;
        args_for_thread[i]->k=k; 
        if(i==0)
        args_for_thread[i]->cnt=(k+1);
        else
        {
          args_for_thread[i]->cnt=k+(args_for_thread[i]->thread_id)*chunk_size;
        }

       }
      for(i = 0; i < NUM_THREADS; i++) 
       pthread_create(&thread_id[i], NULL,division, (void *)args_for_thread[i]);

        for(i = 0; i < NUM_THREADS; i++)
        pthread_join(thread_id[i], NULL);
            
         for(i = 0; i < NUM_THREADS; i++)
         free((void *)args_for_thread[i]);
   */
         for (int j = (k + 1); j < num_elements; j++)
    {     // Reduce the current row
    if (U.elements[num_elements * k + k] == 0)
      {
        printf
    ("Numerical instability detected. The principal diagonal element is zero. \n");
       // return 0;
      }
    U.elements[num_elements * k + j] = (float) (U.elements[num_elements * k + j] / U.elements[num_elements * k + k]);  // Division step
    } 


       U.elements[num_elements * k + k] = 1;
         //  printf("after elimination\n");
          
        for(i = 0; i < NUM_THREADS; i++){
        args_for_thread[i] = (ARGS_FOR_THREAD *)malloc(sizeof(ARGS_FOR_THREAD));
        args_for_thread[i]->thread_id = i; // Provide thread ID  
        args_for_thread[i]->U=&U;
        args_for_thread[i]->k=k; 
        
        if(i==0)
        args_for_thread[i]->cnt=(k+1);
        else
        {
          args_for_thread[i]->cnt=k+(args_for_thread[i]->thread_id)*chunk_size;
        }

      }
         for(i = 0; i < NUM_THREADS; i++) 
       pthread_create(&thread_id[i], NULL,elimination, (void *)args_for_thread[i]);
          
        for(i = 0; i < NUM_THREADS; i++)
        pthread_join(thread_id[i], NULL);
/*
       for(m=0;m<num_elements;m++)
            {
              for(int n=0;n<num_elements;n++)
              printf("   %f   ",U.elements[m*num_elements+n]);
              printf("\n");
            }*/
   //         
         for(i = 0; i < NUM_THREADS; i++)
         free((void *)args_for_thread[i]);
       // printf("after each step\n");*/
        /*printf(" After pthread elimination at %d\n",k);
    for(int m=0;m<num_elements;m++)
            {
              for(int n=0;n<num_elements;n++)
              printf("   %f   ",U.elements[m*num_elements+n]);
              printf("\n");
            }*/   
      }
   /*    printf(" After pthread elimination at %d\n",k);
    for(int m=0;m<num_elements;m++)
            {
              for(int n=0;n<num_elements;n++)
              printf("   %f   ",U.elements[m*num_elements+n]);
              printf("\n");
            }*/
         
}
/*void *division(void *args)
{
    ARGS_FOR_THREAD *args_for_me = (ARGS_FOR_THREAD *)args;
    
    int i,j,num_elements=args_for_me->U->num_rows;
  int chunk_size = (int)floor((float)(num_elements)/(float)NUM_THREADS);
   //int 

//  printf("thread_id %d\n",args_for_me->thread_id );
 // printf("cnt %d\n",args_for_me->cnt );
  //printf("k %d\n",args_for_me->k );

  
   int k= args_for_me->k;
   int last=chunk_size*(args_for_me->thread_id+1)+k;
  //printf("last %d\n", last);
   if(last > num_elements)
    last=num_elements;

//  printf("tid %d\n",args_for_me->thread_id );
   // printf("k %d\n", k);
  //  printf("cnt %d\n", args_for_me->cnt);

    for (j = args_for_me->cnt; j < last; j++)
     {     // Reduce the current row
     if (args_for_me->U->elements[num_elements * (args_for_me->k) + (args_for_me->k)] == 0)
       {
          ("Numerical instability detected. The principal diagonal element is zero. \n");
        }
    //  printf("tid %d\n",args_for_me->thread_id );
     // printf("k %d\n", k);
      //printf("j %d\n", j);
      //printf("index %d\n",num_elements * k + j );
      args_for_me->U->elements[num_elements * k + j] = (float) (args_for_me->U->elements[num_elements * k + j] / args_for_me->U->elements[num_elements * k + k]);  // Division step
      
      }  
}*/
void *elimination(void *args)
{
  ARGS_FOR_THREAD *args_for_me = (ARGS_FOR_THREAD *)args;
   int i,j,num_elements=args_for_me->U->num_rows;
   int chunk_size = (int)floor((float)(num_elements)/(float)NUM_THREADS);
   int k1=args_for_me->cnt;
   int k=args_for_me->k;
   int last=chunk_size*(args_for_me->thread_id+1)+k;
    if(last > num_elements)
    last=num_elements;

    for (i = k1; i < last; i++)
    {
      for (j = k+1 ; j < num_elements; j++)
      args_for_me->U->elements[num_elements * i + j] = args_for_me->U->elements[num_elements * i + j] - (args_for_me->U->elements[num_elements * i + k] * args_for_me->U->elements[num_elements * k + j]);  // Elimination step

      args_for_me->U->elements[num_elements * i + k] = 0;
    }
     /* ARGS_FOR_THREAD *args_for_me = (ARGS_FOR_THREAD *)args;
      int i, j, k, first, last, chunk, extra, tid, num_elements;
      tid = args_for_me->thread_id;
      k = args_for_me->k;
      //U = my_data->upper;
      num_elements = args_for_me->U->num_rows;
      chunk = (num_elements - (k+1))/NUM_THREADS;
      extra = (num_elements - (k+1))%NUM_THREADS;
      first = (tid * chunk) + extra + k;
      last = (tid + 1) * chunk + extra + k;

if(tid == 0)
{
  chunk = chunk + extra;
  first = first - extra;
}

for (i = (first + 1); i <= last; i++)
  {
    for (j = (k + 1); j < num_elements; j++)
      {
        args_for_me->U->elements[num_elements * i + j] =  args_for_me->U->elements[num_elements * i + j] - (args_for_me->U->elements[num_elements * i + k] * args_for_me->U->elements[num_elements * k + j]);  // Elimination step
      
      }
        args_for_me->U->elements[num_elements * i + k] = 0;
  }*/
 
}

/* Function checks if the results generated by the single threaded and multi threaded versions match. */
int
check_results (float *A, float *B, unsigned int size, float tolerance)
{
  for (int i = 0; i < size; i++)
    if (fabsf (A[i] - B[i]) > tolerance)
      return 0;
  return 1;
}


/* Allocate a matrix of dimensions height*width
 * If init == 0, initialize to all zeroes.  
 * If init == 1, perform random initialization. 
*/
Matrix
allocate_matrix (int num_rows, int num_columns, int init)
{
  Matrix M;
  M.num_columns = M.pitch = num_columns;
  M.num_rows = num_rows;
  int size = M.num_rows * M.num_columns;

  M.elements = (float *) malloc (size * sizeof (float));
  for (unsigned int i = 0; i < size; i++)
    {
      if (init == 0)
  M.elements[i] = 0;
      else
  M.elements[i] = get_random_number (MIN_NUMBER, MAX_NUMBER);
    }
  return M;
}


/* Returns a random floating-point number between the specified min and max values. */ 
float
get_random_number (int min, int max)
{
  return (float)
    floor ((double)
     (min + (max - min + 1) * ((float) rand () / (float) RAND_MAX)));
}

/* Performs a simple check on the upper triangular matrix. Checks to see if the principal diagonal elements are 1. */
int
perform_simple_check (const Matrix M)
{
  for (unsigned int i = 0; i < M.num_rows; i++)
    if ((fabs (M.elements[M.num_rows * i + i] - 1.0)) > 0.001)
      return 0;
  return 1;
}