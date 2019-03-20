/*
  This program solves the following equation using ADI method
  u_t = nu*(u_xx + u_yy) + f(u, x, y, t)
  can change f for different reaction function but here we use 0
  INPUT: N int number of grid points
         dt double time step size
         T double final time
         nu double model parameter
   This program applies CUDA.
   AUTHOR: Yuxin Chen
*/

#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <math.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include <curand_kernel.h>
#include "magma.h"
#include "magma_types.h"
#include "magma_lapack.h"

#ifndef M_PI
#define M_PI 3.1415926535897932384626433832795
#endif

const int num_threads = 1024; // N = 1024, num_threads should divide N

__device__ static int dev_N;
__device__ static double dev_dx;
//__device__ static double dev_dy;
__device__ static double dev_nu;
__device__ static double dev_dt;

static void HandleError( cudaError_t err, const char* file, int line) {
    if (err != cudaSuccess) {
        printf("%s in %s at line %d\n", cudaGetErrorString(err),file,line);
        exit(1);
    }
}
#define HANDLE_ERROR( err ) (HandleError(err, __FILE__, __LINE__))


__global__ void IC(double* uold, double* unew)
{
    if (blockIdx.x < dev_N && blockIdx.y < dev_N){
        int ij = blockIdx.x + blockIdx.y * gridDim.x;
        double x = -1 + blockIdx.x*dev_dx;
        double y = -1 + blockIdx.y*dev_dx;
        uold[ij] = -x*y + cos(11*M_PI/2.0*x) * sin(8.0*M_PI*y);
        unew[ij] = -x*y + cos(11*M_PI/2.0*x) * sin(8.0*M_PI*y);
    }
}

// kernel function to transpose the data
__global__ void transpose(double* uold, double* unew){
    unew[threadIdx.x*dev_N + blockIdx.y] = uold[threadIdx.x + blockIdx.y * dev_N];
}

// kernel function to update step A, in x-direction
__global__ void stepAC(double* uold, double* unew)
{
    __shared__ double localu[num_threads];
    
    // global index
    int g_i = (threadIdx.x);// x
    int g_j = (blockIdx.y);// y
    
    int g_ind = g_i + g_j*dev_N; // global 1D index
    
    // load data from global to local memory
    int l_i = threadIdx.x;
    localu[l_i] = uold[g_ind];
    
    __syncthreads();
    
    if (g_i > 0 && g_i < dev_N-1 && g_j > 0 && g_j < dev_N-1) {
        unew[g_ind] = localu[g_i] + dev_nu*dev_dt/dev_dx/dev_dx/2.0*(localu[g_i-1]-2*localu[g_i]+localu[g_i+1]);
    }
}


int main(int argc, char* argv[])
{
    
    // start reading input data using function fscanf here
    int N; // number of grid points in space
    double dt; // time step size
    double T; // terminal time
    double nu;
    
    FILE* inputfile;
    inputfile = fopen(argv[1], "r");
    if (!inputfile) {
        printf("Unable to open input file\n");
        return 0;
    }
    fscanf(inputfile, "%d", &N);
    fscanf(inputfile, "%lf", &dt);
    fscanf(inputfile, "%lf", &T);
    fscanf(inputfile, "%lf", &nu);
    printf("N = %d, dt = %g\n", N, dt);
    fclose(inputfile);
    
    // choose the GPU card
    cudaDeviceProp prop;
    int dev;
    memset(&prop, 0, sizeof(cudaDeviceProp));
    prop.multiProcessorCount = 13;
    HANDLE_ERROR(cudaChooseDevice(&dev, &prop));
    HANDLE_ERROR(cudaSetDevice(dev));
    
    // get maximum thread count for the device
    HANDLE_ERROR(cudaGetDeviceProperties(&prop, dev));
    
    // start recording the time
    cudaEvent_t start, stop;
    HANDLE_ERROR(cudaEventCreate(&start));
    HANDLE_ERROR(cudaEventCreate(&stop));
    HANDLE_ERROR(cudaEventRecord(start, 0));
    
    
    // Addresses into host memory and allocate the memory on the host
    double dx = 2.0/(N-1); // space step size
    double x[N];
    x[0] = -1.0;
    for (int ii = 0; ii < N-1; ii++){
        x[ii+1] = x[ii] +  dx;
    }
    
    // Initiate the MAGMA system
    magma_init();
    
    // temporary data required by dgesv
    magma_int_t *piv, info;
    
    // the right hand side size
    magma_int_t n = N-2;
    magma_int_t m = N;
    
    // error code for MAGMA functions
    magma_int_t err;
    
    // The matrix A and the solution u0 and u1
    double* A;
    double* u0;
    // Addresses into device memory
    double* dev_u0;
    double* dev_u1;
    double* dev_A;
    
    // Allocate matrices on the host using pinned memory
    err = magma_dmalloc_cpu(&A, m*m);
    if (err) {
        printf("oops, an error in memory allocation !\n");
        exit(1);
    }
    err = magma_dmalloc(&dev_A, m*m);
    if (err) {
        printf("oops, an error in memory allocation !\n");
        exit(1);
    }
    err = magma_dmalloc_cpu(&u0, m*m);
    if (err) {
        printf("oops, an error in memory allocation !\n");
        exit(1);
    }
    err = magma_dmalloc(&dev_u0, m*m);
    if (err) {
        printf("oops, an error in memory allocation !\n");
        exit(1);
    }
    err = magma_dmalloc(&dev_u1, m*m);
    if (err) {
        printf("oops, an error in memory allocation !\n");
        exit(1);
    }
    
    // Create temporary storage for the pivots.
    piv = (magma_int_t*)malloc(m*sizeof(magma_int_t));
    
    // make matrix A
    double a = nu*dt/2.0/dx/dx;
    for (int i = 0; i < N; i++){
        for (int j = 0; j < N; j++){
            if (i == j){
                A[i + N*j] = 1 + 2*a;
            }
            else if (i == j+1){
                A[i + N*j] = -a;
            }
            else if (i == j-1){
                A[i + N*j] = -a;
            }
            else {
                A[i + N*j] = 0;
            }
        }
    }
    A[0] = 1;
    A[N*N-1] = 1;
    A[N] = 0;
    A[N*(N-1) - 1] = 0;
    
    // copy the input data from the host to the device
    HANDLE_ERROR(cudaMemcpyToSymbol(dev_dx, &dx, sizeof(double)));
    HANDLE_ERROR(cudaMemcpyToSymbol(dev_dt, &dt, sizeof(double)));
    HANDLE_ERROR(cudaMemcpyToSymbol(dev_nu, &nu, sizeof(double)));
    HANDLE_ERROR(cudaMemcpyToSymbol(dev_N, &N, sizeof(int)));
    
    // Initialize the function data, use 2 dim block
    dim3 initDim(N,N);
    IC<<<initDim, 1>>>(dev_u0, dev_u1);
    
    // copy initial data back to host
    magma_dgetmatrix(m, m, dev_u0, m, u0, m);
    
    // open the file to write output
    FILE* outputfile;
    outputfile = fopen(argv[2], "wb");
    fwrite(&N, sizeof(int), 1, outputfile);
    fwrite(x, sizeof(double), N, outputfile);
    fwrite(u0, sizeof(double), N*N, outputfile); // write out the initial condition
    
    // using num_blocks blocks and each block contains num_threads*num_threads many threads
    dim3 meshDim(1,N);
    dim3 thrddim(num_threads,1);
    
    // copy matrix on host onto device
    magma_dsetmatrix(m, m, A, m, dev_A, m);
    magma_dgetrf_gpu(m, m, dev_A, m, piv, &info);
    
    // time loop starts
    int iter = round(T/dt);
    int piter = round(T/dt/4.0);
    for (int k = 1; k <= iter; k++){
        
        // step A
        stepAC<<<meshDim, thrddim>>>(dev_u0, dev_u1);
        transpose<<<meshDim, thrddim>>>(dev_u1, dev_u0);
        
        // step B, call implicit solver
        magma_dgetrs_gpu(MagmaNoTrans, m, n, dev_A, m, piv, &(dev_u0[N]), m, &info);
        
        // step C
        stepAC<<<meshDim, thrddim>>>(dev_u0, dev_u1);
        
        // transpose data
        transpose<<<meshDim, thrddim>>>(dev_u1, dev_u0);
        
        // step D, call implicit solver
        magma_dgetrs_gpu(MagmaNoTrans, m, n, dev_A, m, piv, &(dev_u0[N]), m, &info);
 
        // write data to outputfile
        if (k%piter == 0){
            magma_dgetmatrix(m, m, dev_u0, m, u0, m);
            fwrite(u0, sizeof(double), N*N, outputfile);
        }
    }
    
    // close outputfile
    fclose(outputfile);
    
    // clean up all the allocated memory
    free(u0);
    free(A);
    free(piv);
    magma_free(dev_u0);
    magma_free(dev_u1);
    magma_free(dev_A);
    
    HANDLE_ERROR(cudaEventRecord(stop, 0));
    HANDLE_ERROR(cudaEventSynchronize(stop));
    float elapsedTime;
    HANDLE_ERROR(cudaEventElapsedTime(&elapsedTime, start, stop));
    
    HANDLE_ERROR(cudaEventDestroy(start));
    HANDLE_ERROR(cudaEventDestroy(stop));
    
    printf("Execution time = %e ms\n", elapsedTime);
    magma_finalize();
    return 0;
}

