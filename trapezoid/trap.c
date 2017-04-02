/*  Purpose: Calculate definite integral using trapezoidal rule.
 *
 * Input:   a, b, n
 * Output:  Estimate of integral from a to b of f(x)
 *          using n trapezoids.
 *
 * Compile: gcc -o trap trap.c -lpthread -lm
 * Usage:   ./trap
 *
 * Note:    The function f(x) is hardwired.
 *
 */

#ifdef _WIN32
#  define NOMINMAX 
#endif

// includes, system
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <time.h>
#include <sys/time.h> 
#include <pthread.h>

#define LEFT_ENDPOINT 5
#define RIGHT_ENDPOINT 1000
#define NUM_TRAPEZOIDS 100000000
#define NUM_THREADS 4

double compute_using_pthreads(float, float, int, float);
double compute_gold(float, float, int, float);
void *	partial_integral_func(void *);
typedef struct args_for_thread_s{
        int thread_id; // The thread ID
		float a;
		float h;
		int chunk_size;
		double *partial_integral; 
	//	int total_num;
} ARGS_FOR_THREAD;

int main(void) 
{
	int n = NUM_TRAPEZOIDS;
	float ts,tp;
	float a = LEFT_ENDPOINT;
	float b = RIGHT_ENDPOINT;
	float h = (b-a)/(float)n; // Height of each trapezoid  
	printf("The height of the trapezoid is %f \n", h);
	 struct timeval start, stop;
     gettimeofday (&start, NULL);
 	double reference = compute_gold(a, b, n, h);
	gettimeofday (&stop, NULL);
  printf ("CPU run time = %0.2f s. \n",
    (float) (stop.tv_sec - start.tv_sec +
       (stop.tv_usec - start.tv_usec) / (float) 1000000));
ts=(float) (stop.tv_sec - start.tv_sec +
       (stop.tv_usec - start.tv_usec) / (float) 1000000);
     
   printf("Reference solution computed on the CPU = %f \n", reference);

	/* Write this function to complete the trapezoidal on the GPU. */
       gettimeofday (&start, NULL);
	double pthread_result = compute_using_pthreads(a, b, n, h);
	gettimeofday (&stop, NULL);
  printf ("Pthread implementation run time = %0.2f s. \n",
    (float) (stop.tv_sec - start.tv_sec +
       (stop.tv_usec - start.tv_usec) / (float) 1000000));
tp=(float) (stop.tv_sec - start.tv_sec +
       (stop.tv_usec - start.tv_usec) / (float) 1000000);
    
	printf("Solution computed using pthreads = %f \n", pthread_result);
	printf("Speedup: %f\n",ts/tp );

} 


/*------------------------------------------------------------------
 * Function:    f
 * Purpose:     Compute value of function to be integrated
 * Input args:  x
 * Output: (x+1)/sqrt(x*x + x + 1)

 */
float f(float x) {
		  return (x + 1)/sqrt(x*x + x + 1);
}  /* f */

/*------------------------------------------------------------------
 * Function:    Trap
 * Purpose:     Estimate integral from a to b of f using trap rule and
 *              n trapezoids
 * Input args:  a, b, n, h
 * Return val:  Estimate of the integral 
 */
double compute_gold(float a, float b, int n, float h) {
   double integral;
   int k;
   integral = (f(a) + f(b))/2.0;
   for (k = 1; k <= n-1; k++) {
     integral += f(a+k*h);
   }
   integral = integral*h;
   return integral;
}  

/* Complete this function to perform the trapezoidal rule on the GPU. */
double compute_using_pthreads(float a, float b, int n, float h)
{
		 pthread_t thread_id[NUM_THREADS]; // Data structure to store the thread IDs
    	pthread_attr_t attributes; // Thread attributes
    	pthread_attr_init(&attributes); // Initialize the thread attributes to the default values 
		int i;
		ARGS_FOR_THREAD *args_for_thread[NUM_THREADS];
		double *partial_integral = (double *)malloc(sizeof(double) * NUM_THREADS);
		for(i=0;i<NUM_THREADS;i++)
		{
			args_for_thread[i] = (ARGS_FOR_THREAD *)malloc(sizeof(ARGS_FOR_THREAD));
			args_for_thread[i]->thread_id=i;
			args_for_thread[i]->a=a;
			args_for_thread[i]->h=h;
			args_for_thread[i]->chunk_size=(int)floor((float)(n)/(float)NUM_THREADS);
			args_for_thread[i]->partial_integral=partial_integral;
			//args_for_thread[i]->total_num=NUM_TRAPEZOIDS;
		}		
		for(i = 0; i < NUM_THREADS; i++)
        pthread_create(&thread_id[i], &attributes, partial_integral_func, (void *)args_for_thread[i]);
					 
    /* Wait for the workers to finish. */
    for(i = 0; i < NUM_THREADS; i++)
        pthread_join(thread_id[i], NULL);
		 
    /* Accumulate partial sums. */ 
    double integral_final =0.0 ; 

    for(i = 0; i < NUM_THREADS; i++)
        integral_final += partial_integral[i];

	integral_final=integral_final+(f(a)+f(b))/2.0;	
    /* Free dynamically allocated data structures. */
    free((void *)partial_integral);
    for(i = 0; i < NUM_THREADS; i++)
        free((void *)args_for_thread[i]);	

		  return (integral_final*h);
}
void *partial_integral_func(void* args)
{
ARGS_FOR_THREAD *args_for_me= (ARGS_FOR_THREAD *)args;
int k;
float a=args_for_me->a;
float h=args_for_me->h;
int chunk_size=args_for_me->chunk_size;

int first=(args_for_me->thread_id*chunk_size)+1;
printf("Tid %d\n",args_for_me->thread_id );
int last=first+chunk_size;

double  integral_me=0.0;
float temp=0.0;
//if(last>args_for_me->total_num)
//	last=args_for_me->total_num;

for (k = first; k <= last; k++) {
    temp=(a+k*h);
    integral_me =  integral_me + f(temp);
   }
  // printf("partial sum  %f \n",integral_me );
args_for_me->partial_integral[args_for_me->thread_id]=integral_me;

}


