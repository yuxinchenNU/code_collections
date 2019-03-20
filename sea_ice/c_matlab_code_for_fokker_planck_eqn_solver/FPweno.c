/*
    Solve Fokker-Planck equation with semi-implicit method
    implicit in diffusion (Crank-Nicolson) (2-2) scheme
    explicit in advection-reaction terms, weno reconstructed stencils
    applied the method in 
    "High-order Finite Difference and Finite Volume WENO Schemes and Discontinuous Galerkin Methods for CFD"
    to solve fokker-planck equation that describes probability density function of
    random walk with a time-varying external force.
    Used in the sea ice project
 
*/

#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <math.h>
#include <Accelerate/Accelerate.h>

// declare function used in the model
double drift(double,double,double,double);

int main(int argc, char* argv[])
{
    // read input parameters
	if (argc != 3) {
		printf("Incorrect usage: only enter the input data file name\n");
		return 0;
	}   
	FILE* inputfile = fopen(argv[1], "r");
	if (!inputfile) {
		printf("Unable to open input file\n");
		return 0;
	}
	// start reading input data using function fscanf here
	int N; // number of grid points in space
	double dt; // time step size
	double T; // terminal time
    double ic; // starting point
    double L; // x\in [-L L]
    double A; // amplitude of periodic forcing
    double eps; //
    double sigma; // noise intensity

	fscanf(inputfile, "%d", &N);
	fscanf(inputfile, "%lf", &dt); 
	fscanf(inputfile, "%lf", &T);
    fscanf(inputfile,"%lf", &ic);// initial position
    fscanf(inputfile,"%lf", &L); // x\in [-L, L]
    fscanf(inputfile,"%lf", &A);
    fscanf(inputfile,"%lf", &eps);
    fscanf(inputfile,"%lf", &sigma);

	fclose(inputfile); // now we have all the parameters we want, close inputfile

    int N0 = N; // keep the original number of grid points in space
    N = N0 + 2; // add two ghost point
    
	clock_t start = clock(); // clock starts ticking

	// additional parameters
	int M = round(T/dt); // number of grid points in time
	double D = sigma*sigma/2;
	double dx = 2.0*L/(N0-1); // space step size
    double constant = D*dt/2.0/dx/dx; // simplifed symbol used in CN step
    
	// define space vector x and time vector t
	double x[N]; // including two ghost points
    double t[M];
	x[0] = -L-dx; // since this is the ghost point to the leftv
	for (int ii = 0; ii < N-1; ii++){
		x[ii+1] = x[ii] +  dx;
	}
    t[0] = 0;
    for (int ii = 0; ii < M; ii++){
        t[ii+1] = t[ii] +  dt;
    }
    
    // define face vector xh
    double xh[N];
    xh[0] = x[0] + 0.5*dx;
    for (int ii = 0; ii < N-1; ii++){
        xh[ii+1] = xh[ii] + dx;
    }
	

	// allocate memory
    double (*p)[N] = malloc(M*sizeof(*p));
    double u[N];
    // initialization step
    double s = 0.1; // to approximate delta function
    for (int space = 0; space < N; space++)
    {
        p[0][space] = 1/sqrt(2*M_PI)/s*exp(-pow(x[space]-ic,2)/2/s/s);
        u[space] = p[0][space];
    }
    
    
    // parameters passed in banded matrix solver in LAPACK
    int info;
    double* b = malloc(N*sizeof(double)); // right hand side of the linear system

    double (*Mat)[N] = malloc(N*sizeof(*Mat));
       int piv[N];
    int one = 1;
    
    // save a copy of the original data
//    double* borig = malloc(N*sizeof(double));
//    double (*MatOrig)[N] = malloc(N*sizeof(*MatOrig));
    
    // fill in the matrix
    for (int i = 1; i < N-1; i++){
        Mat[i][i] = 1+2*constant;
        Mat[i+1][i] = -constant;
        Mat[i-1][i] = -constant;
    }
    Mat[0][0] = D/2.0/dx;
    Mat[2][0] = -D/2.0/dx;
    Mat[N-1][N-1] = -D/2.0/dx;
    Mat[N-3][N-1] = D/2.0/dx;
 
    
    // allocate memories for some variables used in WENO construction
    double vm[N]; // fhat_minus at x_{i+1/2}
    double vp[N]; // fhat_plus at x_{i+1/2}
    double fhat[N];
    double fu[N];
    double fuAppr1[N];
    double fuAppr2[N];
    double fuAppr3[N];
    double fuApprt1[N];
    double fuApprt2[N];
    double fuApprt3[N];
    double a[N]; // Roe speed
    
    // some constants used in WENO construction
    double c[9] = {1.0/3.0, 5.0/6.0, -1.0/6.0, -1.0/6.0, 5.0/6.0, 1.0/3.0, 1.0/3.0, -7.0/6.0, 11.0/6.0};
    double ct[9] = {11.0/6.0, -7.0/6.0, 1.0/3.0, 1.0/3.0, 5.0/6.0, -1.0/6.0, -1.0/6.0, 5.0/6.0, 1.0/3.0};
    double c0[6] = {1.5, -0.5, 0.5, 0.5, -0.5, 1.5};
    double d1 = 2.0/3.0;
    double d2 = 1.0/3.0; // used for computing vm(2)
    double v1, v2, beta1, beta2, alpha1, alpha2, w1, w2; // for computing vm(2)
    double v1t, v2t, alphat1, alphat2, wt1, wt2; // for computing vm(2)
    double dt1 = 1.0/3.0; // used for computing vp(N-2)
    double dt2 = 2.0/3.0;
    double d[3] = {3.0/10.0, 3.0/5.0, 1.0/10.0};
    double dT[3] = {1.0/10.0, 3.0/5.0, 3.0/10.0};
    double beta[3];
    double alpha[3];
    double w[3];
    double wt[3];
    double alphat[3];
    double epsilon = 0.000001;
    
    // model steps
    for (int n = 1; n <= M; n++){
        
        //////// advection term update
        
        // vector to store f(u_i), used to compute h^bar_i
        for (int i = 0; i < N; i++){
            fu[i] = drift(x[i],t[n],eps,A)*u[n];
        }
        /// step1, compute the approximation from different neighboring stencils
        for (int i = 2; i < N-2; i++){
            // v(r)'s eqn 2.51
            fuAppr3[i] = c[0]*fu[i-2] + c[1]*fu[i-1] + c[2]*fu[i];
            fuAppr2[i] = c[3]*fu[i-1] + c[4]*fu[i] + c[5]*fu[i+1];
            fuAppr1[i] = c[6]*fu[i] + c[7]*fu[i+1] + c[8]*fu[i+2];
        }
        
        // step2, obtain the k reconstructed values v_{i-1/2}^(r) using (2.10) based on stencils (2.50)
        for (int i = 2; i < N-2; i++){
            // v(r)_{i-1/2}
            // fuApprt only has information at xh(2:N-3)
            fuApprt1[i-1] = ct[0]*fu[i] + ct[1]*fu[i+1] + ct[2]*fu[i+2];
            fuApprt2[i-1] = ct[3]*fu[i-1] + ct[4]*fu[i] + ct[5]*fu[i+1];
            fuApprt3[i-1] = ct[6]*fu[i-2] + ct[7]*fu[i-1] + ct[8]*fu[i];
        }
        
        // step 3: find nonlinear weights (use the weights given in Jiang & Shu, JCP 1996)
        for (int i = 2; i < N-2; i++){
            beta[2] = 13.0/12.0*pow((fu[i-2] - 2*fu[i-1] + fu[i]), 2) + 0.25*pow((fu[i-2] - 4*fu[i-1] + 3*fu[i]), 2); // 2.63
            beta[1] = 13.0/12.0*pow((fu[i-1] - 2*fu[i] + fu[i+1]), 2) + 1.0/4*pow((fu[i-1] - fu[i+1]), 2);
            beta[0] = 13.0/12*pow((fu[i] - 2*fu[i+1] + fu[i+2]), 2) + 1.0/4.0*pow((3*fu[i] - 4*fu[i+1] + fu[i+2]), 2);
            alpha[0] = d[0]/pow((epsilon + beta[0]), 2); // eqn 2.59
            alpha[1] = d[1]/pow((epsilon + beta[1]), 2);
            alpha[2] = d[2]/pow((epsilon + beta[2]), 2);
            w[0] = alpha[0]/(alpha[0]+alpha[1]+alpha[2]); // eqn 2.58
            w[1] = alpha[1]/(alpha[0]+alpha[1]+alpha[2]);
            w[2] = alpha[2]/(alpha[0]+alpha[1]+alpha[2]);
            // step 4 in procedure 2.2
            alphat[0] = dT[0]/pow((epsilon + beta[0]), 2);
            alphat[1] = dT[1]/pow((epsilon + beta[1]), 2);
            alphat[2] = dT[2]/pow((epsilon + beta[2]), 2);
            wt[0] = alphat[0]/(alphat[0]+alphat[1]+alphat[2]);
            wt[1] = alphat[1]/(alphat[0]+alphat[1]+alphat[2]);
            wt[2] = alphat[2]/(alphat[0]+alphat[1]+alphat[2]);
            // compute vminus{i+1/2} step 5
            // only has information at xh(3:N-2), need 2 stencils to compute
            // vm(2)
            vm[i] = w[0]*fuAppr1[i] + w[1]*fuAppr2[i] + w[2]*fuAppr3[i];
            // compute vplus{i-1/2} step 5
            // only has information at xh(2:N-3), need 2 stencils to compute
            // vp(N-2)
            vp[i-1] = wt[0]*fuApprt1[i-1] + wt[1]*fuApprt2[i-1] + wt[2]*fuApprt3[i-1];
        }
        
        ////// compute vm(2), vm[1] in C
        // compute v(r)_{i+1/2}, i = 2
        v1 = c0[2]*fu[1] + c0[3]*fu[2]; // r = 1;
        v2 = c0[4]*fu[0] + c0[5]*fu[3]; // r = 2;
        // compute vm(2), vm[1] in C
        beta1 = pow((fu[2] - fu[1]), 2);
        beta2 = pow((fu[1] - fu[0]), 2);
        alpha1 = d1/pow((epsilon + beta1),2);
        alpha2 = d2/pow((epsilon + beta2),2);
        w1 = alpha1/(alpha1 + alpha2);
        w2 = alpha2/(alpha1 + alpha2);
        vm[1] = w1*v1 + w2*v2;
        
        
        /// compute vp(N-2), vp[N-3]
        v1t = c0[0]*fu[N-2] + c0[2]*fu[N-1];
        v2t = c0[2]*fu[N-3] + c0[3]*fu[N-2];
        alphat1 = dt1/pow((epsilon + beta1),2);
        alphat2 = dt2/pow((epsilon + beta2),2);
        wt1 = alphat1/(alphat1 + alphat2);
        wt2 = alphat2/(alphat1 + alphat2);
        vp[N-3] = wt1*v1t + wt2*v2t;
        
        /// compute vm(1) (vm[0] in C) at xh(1) (xh[0] in C), ie x(1+1/2) where x(1) is the ghost point
        vm[0] = fu[0]; // w1 = d1 = 1, c11 = 1;
        vp[0] = fu[1];
        
        /// compute vm(N-1) (vm[N-2] in C) and vp(N-1) (vp[N-2] in C) at xh(N-1), ie x(N-1+1/2) where x(N) is the ghost point
        vm[N-2] = fu[N-2]; // w1 = d1 = 1, c11 = 1
        vp[N-2] = fu[N-1];
        
        /// compute Roe speed
        for (int i = 0; i < N-1; i++){
            a[i] = (fu[i+1] - fu[i])/(u[i+1] - u[i]);
        }
        
        /// loop through the faces starting from 1 to N-1, (0 to N-2 in C)
        for (int i = 0; i < N-2; i++) {
            if (a[i] > 0)
                fhat[i] = vm[i];
            else
                fhat[i] = vp[i];
        }
        
        /////////// DIFFUSION TERM UPDATE  /////////////
        // continue filling in the matrix
        Mat[1][0] = drift(x[1],t[n],eps,A);
        Mat[N-2][N-1] = drift(x[N-2],t[n],eps,A);
        // update right hand side
        for(int i = 1; i < N; i++){
            b[i] = u[i] - dt/dx*(fhat[i] - fhat[i-1]) + constant*(u[i+1] - 2*u[i] + u[i-1]);
        }
        dgesv_(&N, &one, &(Mat[0][0]), &N, piv, b, &N, &info);
        
        for (int i = 0; i < N; i++){
            u[i] = b[i];
        }
    }
    
	
	printf("%d\n",N);
	printf("Time elapsed: %g seconds\n", (float)(clock()-start)/CLOCKS_PER_SEC);
	// read file 
	FILE* outputfile = fopen(argv[2], "w");
	
    fwrite(&N,sizeof(int),1,outputfile);
    fwrite(&M,sizeof(int),1,outputfile);
    fwrite(x,sizeof(double),N,outputfile);
    fwrite(u, sizeof(double), N, outputfile);
    fclose(outputfile);
	free(p);
    free(b);
    free(Mat);

	return 0;



}

double drift(double x, double t, double eps, double A)
{
	double v = x - pow(x,3) + A*cos(eps*2.0*M_PI*t);
	return v;
}






