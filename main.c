#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <mpi.h>
#include <petsc.h>
#include <petscksp.h>
#include <petscvec.h>
#include <petscmat.h>
#include <petscsys.h>
#include "main.h"
#include "memory.h"
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define RECOVER
#define MOVE
#define TEMP
#define WRITE

int main(int argc, char *argv[]){

	PetscInitialize(&argc, &argv, 0, 0);
	int rank, nbproc;
	MPI_Comm_rank(PETSC_COMM_WORLD, &rank );
	MPI_Comm_size(PETSC_COMM_WORLD, &nbproc);
	printf("Hello from rank %d \n", rank);

	/* DIMENSIONS */
	Dp = 1;
	ratio_L_d = 1;
	ratio_d_Dp = 30;
	ratio_Dp_h = 30;
	d = ratio_d_Dp*Dp;
	H = d/2.; /* half -height of the channel */
	L = ratio_L_d*d; /* length of the channel is 'ratio' times its width */

	/* FLOW */
	Rp = 40.;
	Um = Rp*nu/Dp; /* By definition of Rp */

	Tm0 = 100. + 273.15; // cup-mixing temperature at the inlet
	Tp0 = 20. + 273.15;

	Ns = 2; /* Number of species */
	Np = 1; /* Number of particles */

	/* PHYSICAL PARAMETERS */
	rho_r = rho_p/rho_f;
	cr = cp/cf;
	alpha_f = nu/Pr;

	/* SPECIES PARAMETERS */
	double CA0 = 1.;
	double CB0 = 0.;
	Df = make1DDoubleArray(Ns);
	Df[0] = 2e-5;
	Df[1] = 2e-5;

	/* GRID */
	N = ratio_Dp_h*ratio_d_Dp;
	M = ratio_L_d*N;
	n = N + 2; /*for ghost points */
	m = M + 2; /*for ghost points */
	h = d/N;

	/* TIME INTEGRATION */

	double dt_CFL = CFL*h/Um;
	double dt_diff = r*h*h/nu;
	double refine_dt = 10.;
	dt = fmin(dt_CFL, dt_diff)/refine_dt;
	double ratio_dtau_dt = 1e-2;
	dtau = ratio_dtau_dt*dt;

	if(rank == 0){
		printf("Rep = %f\n", Rp);
		printf("ratio L/d = %d \n", ratio_L_d);
		printf("ratio d/dp = %d \n", ratio_d_Dp);
		//printf("Umax = %f\n", Umax);
		printf("Um = %f\n", Um);
		printf("dt_CFL = %f\n", dt_CFL);
		printf("dt_diff = %f\n", dt_diff);
		printf("dt is reduced by a factor %f with respect to CFL \n", refine_dt);
		printf("dt = %f\n", dt);
		printf("dtau = %f\n", dtau);
		printf("ratio dtau/dt = %f \n", ratio_dtau_dt);
	}
	int T_write; 
	int N_write;
	int num;
	if (argc >= 3){
		sscanf(argv[1],"%d",&T_write);
		sscanf(argv[2],"%d",&N_write);
	}
	else{
		T_write = 20; /* number of time steps between two writings */ 
		N_write = 200; /* number of times we write in files */
	}
	double Tf = N_write*T_write*dt;
	int nKmax = 2;
	int Kmax = 20; /* number of ramping steps */

	double t, t_start;
	int iter, iter_start;
	int l, l_start;

	double ramp = 1./Kmax;
	if(rank == 0){
		printf("Write every %d * dt \n", T_write);
		printf("Write %d times \n", N_write);
		printf("Final time : %f \n \n", Tf);
		printf("ramp = %f \n", ramp);
	}

	/* PARTICLES POSITION */
	double* xg = calloc(Np,sizeof(double));
	double* yg = calloc(Np,sizeof(double));
	double* dp = calloc(Np,sizeof(double));
	double* rp = calloc(Np,sizeof(double));
	double* theta = calloc(Np,sizeof(double));
	double* Sp = calloc(Np,sizeof(double));
	double* II = calloc(Np,sizeof(double));

	xg[0]=H;
	yg[0]=H;
	dp[0]=Dp;
	rp[0]=dp[0]/2.;
	theta[0] = 0;

	for(int k=0; k<Np; k++){
		Sp[k]=M_PI*rp[k]*rp[k];
		II[k]=(dp[k]*dp[k])/8.; /* IN 2-D !! */
	}

#ifdef WRITE
	DEFINE_FILES
		if(rank==0)
		{
			/* Open files to write the results */
			OPEN_FILES	
			/* Write DATA */
			WRITE_DATA
		}
#endif

	/* Allocate memory */
	MEMORY_ALLOCATION

	/* Try to recover previous state */
	OPEN_STATE

#ifndef RECOVER
		state_file = NULL;
#endif
	if(state_file){
		RECOVER_STATE
	}
	else{

		/**------------------------------- INITIALIZATION of the domain -------------------------------**/

		/* VELOCITY : horizontal flow Um  */
		for(int i=0; i<m; i++){
			for(int j=0; j<n; j++){
				U[i][j] = Um; /*Plug-Flow*/
				Ustar[i][j] = U[i][j];
				Told[i][j] = Tm0;
			}
		}

		/*Initalization of particles temperatures */
		for(int k=0; k<Np; k++){
			Tp[k] = Tp0;
		}

		/* ------------------------- FIRST STEP WITH EULER EXPLICIT SCHEME ---------------------------  */
		/* Location PARTICLES */
		/* GET MASK FUNCTIONS */
		get_masks(Ip_S, Ip_U, Ip_V, I_S, I_U, I_V, xg, yg, rp, theta, coloring);

		/* Initialize solid temperature and field based on initial particle temperature/species concentration */
		get_Ts(Ts, Ip_S, Tp);

		/* Get Ghost points */
		get_ghosts(U, V, P, Told, Cold, 0, 0);

		/* Solve the flow */

		get_Ustar_EE(U, Ustar, V, Aold, Us, P, I_U, ramp);
		get_Vstar_EE(V, Vstar, U, Bold, Vs, P, I_V, ramp);

		//old_poisson_solver(Ustar, Vstar, phi, R);
		poisson_solver(Ustar, Vstar, phi, rank, nbproc);
		update_flow(U, V, P, Ustar, Vstar, phi);

		/* Solve temperature and specices */
#ifdef TEMP
		update_temp_species_EE(U, V, T, Told, Ts, C_Told, C, Cold, Cs, C_Cold, I_S, ramp);
#endif

		/** ------------------------------- RAMPING ------------------------------- **/

		/* Starting from a uniform flow without the particles (i.e. ramp = 0 everywhere) and
		   then progressively increasing ramp until ramp*chi = 1 at the position of the ellipse.
		   At the end of the loop we obtain the steady solution of the flow for a given
		   fixed position of the ellipse. */

		for(int K=2; K<=nKmax*Kmax; K++){

			/* role de nKmax : atteindre une solution stable avec les particules fixées avant de les "lacher" */
			ramp = fmin(1.,(double) K/Kmax);

			PetscPrintf(PETSC_COMM_WORLD, "ramp = %f\n",ramp);

			get_ghosts(U, V, P, Told, Cold, 0, 0);

			get_Ustar(U, Ustar, V, A, Aold, Us, P, I_U, ramp);
			get_Vstar(V, Vstar, U, B, Bold, Vs, P, I_V, ramp);

			poisson_solver(Ustar, Vstar, phi, rank, nbproc);
			//old_poisson_solver(Ustar, Vstar, phi, R);
			update_flow(U, V, P, Ustar, Vstar, phi);
#ifdef TEMP
			update_temp_species(U, V, T, Told, Ts, C_T, C_Told, C, Cold, Cs, C_C, C_Cold, I_S, ramp);
#endif

			for (int k=0; k<Np; k++){
				compute_force_torque_fluxes(F, G, M, Qp, Qmp, Ip_U, Ip_V, Ip_S, U, V, T, C, Us, Vs, Ts, Cs, xg, yg, rp, Sp, II, F_drag, F_lift, Torque, Q_heat, Phi_species, k);
			}

		}

#ifdef WRITE
		/*INITIAL SOLUTION AFTER RAMPING */
		if(rank==0){
			fprintf(fichier_position, "%3.13e \t %3.13e \t %3.13e \n ", xg[0], yg[0], theta[0]);
			fflush(fichier_position);
			fprintf(fichier_forces, "%3.13e \t %3.13e \t %3.13e \n ", F_drag[0], F_lift[0], Torque[0]);
			fflush(fichier_forces);
			fprintf(fichier_fluxes, "%3.13e \t %3.13e \n ", Q_heat[0], Phi_species[0]);
			fflush(fichier_fluxes);
			writeFile(fichier_U,U,0,m-1,1,n-1);
			writeFile(fichier_V,V,1,m-1,0,n-1);
			writeFile(fichier_P,P,1,m-1,1,n-1);
			writeFile(fichier_T,T,1,m-1,1,n-1);
			writeFile(fichier_CA,C[0],1,m-1,1,n-1);
			writeFile(fichier_CB,C[1],1,m-1,1,n-1);
			writeFile(fichier_mask,coloring,1,m-1,1,n-1);
			fprintf(fichier_Tp, "%3.13e \n",Tp[0]);
			fflush(fichier_Tp);
		}
#endif
		/** -------------------------------TIME STEPPING FROM BEGINNING ------------------------------- **/
		iter_start = 1;
		l_start = 1;
		t_start = 0.;
	}

	/** -------------------------------TIME STEPPING ------------------------------- **/
	iter = iter_start;
	t = t_start;
	if(rank==0){
		printf("iter %d : t = %f\n", iter, t);
	}
	fflush(stdout);
	while(t < Tf){
		for (int k = 0; k<Np; k++){
			/*Update particles positions and velocities  */
#ifdef MOVE
			update_Xp(xg, yg, theta, Up, Vp, Wp, k);
			update_Up(Up, Vp, Wp, F, G, M, k);
#endif
#ifdef TEMP
			update_Tp(Tp, Qp, k);
			update_Cp(Cp, Qmp, k);
#endif
		}
		/*Compute the mask functions */
		get_masks(Ip_S, Ip_U, Ip_V, I_S, I_U, I_V, xg, yg, rp, theta, coloring);

		/* Deduce solid velocity field */
		get_Us_Vs(Us, Vs, Ip_U, Ip_V, Up, Vp, Wp, xg, yg);

#ifdef TEMP
		get_Ts(Ts, Ip_S, Tp);
		get_Cs(Cs, Ip_S, Cp);
#endif
		get_ghosts(U, V, P, Told, Cold, CA0, CB0);

		get_Ustar(U, Ustar, V, A, Aold, Us, P, I_U, ramp);
		get_Vstar(V, Vstar, U, B, Bold, Vs, P, I_V, ramp);

		clock_t t_init = clock();
		//old_poisson_solver(Ustar, Vstar, phi, R);
		poisson_solver(Ustar, Vstar, phi, rank, nbproc);
		clock_t t_final = clock();
		double t_Poisson = ((double) (t_final - t_init))/CLOCKS_PER_SEC;
		PetscPrintf(PETSC_COMM_WORLD, "Poisson solver took %f seconds \n", t_Poisson);

		update_flow(U, V, P, Ustar, Vstar, phi);
#ifdef TEMP
		update_temp_species(U, V, T, Told, Ts, C_T, C_Told, C, Cold, Cs, C_C, C_Cold, I_S, ramp);
#endif
		int flag_out = 0;  
		/*For each particle :*/
		for (int k = 0; k<Np; k++){
			flag_out += compute_force_torque_fluxes(F, G, M, Qp, Qmp, Ip_U, Ip_V, Ip_S, U, V, T, C, Us, Vs, Ts, Cs, xg, yg, rp, Sp, II, F_drag, F_lift, Torque, Q_heat, Phi_species, k);
		}

		if(flag_out > 0){
			break; 
		}

		/* Increment time */
		t += dt;
		iter++;

		PetscPrintf(PETSC_COMM_WORLD,"\n \n iter %d : t = %f\n", iter, t);
		fflush(stdout);
#ifdef WRITE
		if(iter % T_write == 0 && rank == 0){
			fprintf(fichier_position, "%3.13e \t %3.13e \t %3.13e \n ", xg[0], yg[0], theta[0]);
			fflush(fichier_position);
			fprintf(fichier_forces, "%3.13e \t %3.13e \t %3.13e \n ", F_drag[0], F_lift[0], Torque[0]);
			fflush(fichier_forces);
			fprintf(fichier_fluxes, "%3.13e \t %3.13e \n ", Q_heat[0], Phi_species[0]);
			fflush(fichier_fluxes);
			writeFile(fichier_U,U,0,m-1,1,n-1);
			writeFile(fichier_V,V,1,m-1,0,n-1);
			writeFile(fichier_P,P,1,m-1,1,n-1);
			writeFile(fichier_T,T,1,m-1,1,n-1);
			writeFile(fichier_CA,C[0],1,m-1,1,n-1);
			writeFile(fichier_CB,C[1],1,m-1,1,n-1);
			fprintf(fichier_Tp, "%3.13e \n",Tp[0]);
			fflush(fichier_Tp);
#ifdef MOVE
			writeFile(fichier_mask,coloring,1,m-1,1,n-1);
#endif
		}
#endif
		t += dt;
		iter++;

		PetscPrintf(PETSC_COMM_WORLD,"\n \n iter %d : t = %f\n", iter, t);            
		fflush(stdout);  

	}
	printf("rank %d : job done", rank);

#ifdef WRITE
	/*Save current state */
	SAVE_STATE

		if (rank ==0){
			/* Close files */
			CLOSE_FILES
		}
#endif

	/* Free memory */
	FREE_MEMORY

		PetscFinalize();
	return 0;
}

int compute_force_torque_fluxes(double** dudt, double** dvdt, double** domegadt, double** dTdt, double*** dCdt, double*** Ip_U, double*** Ip_V, double*** Ip_S, double** U, double** V, double** T, double*** C, double** Us, double** Vs,double** Ts, double*** Cs, double* xg, double* yg, double* rp, double* Sp, double* II, double* F_drag, double* F_lift, double* Torque, double* Q_heat, double* Phi_species, int k)
{
	/* Force along x-direction */
	dudt[k][0] = dudt[k][1]; /* n-2*/
	dudt[k][1] = dudt[k][2]; /* n-1*/
	dudt[k][2] = 0.; /* n*/

	/* Force along y-direction */
	dvdt[k][0] = dvdt[k][1];
	dvdt[k][1] = dvdt[k][2];
	dvdt[k][2] = 0.;

	/* Moment along z-direction */\
		domegadt[k][0] = domegadt[k][1];
	domegadt[k][1] = domegadt[k][2];
	domegadt[k][2] = 0.;

	/* Particle heat balance  */
	dTdt[k][0] = dTdt[k][1]; /* n-2*/
	dTdt[k][1] = dTdt[k][2]; /* n-1*/
	dTdt[k][2] = 0.; /* n*/

	for(int s=0; s<Ns; s++){
		dCdt[k][s][0] = dCdt[k][s][1];
		dCdt[k][s][1] = dCdt[k][s][2];
		dCdt[k][s][2] = 0.;
	}

	int startX = floor((xg[k]-rp[k])/h);
	PetscPrintf(PETSC_COMM_WORLD, "startX = %d \t", startX);
	int endX = ceil((xg[k]+rp[k])/h);
	PetscPrintf(PETSC_COMM_WORLD,"endX = %d \t", endX);

	int startY = floor((yg[k]-rp[k])/h);
	PetscPrintf(PETSC_COMM_WORLD,"startY = %d \t", startY);
	int endY = ceil((yg[k]+rp[k])/h);
	PetscPrintf(PETSC_COMM_WORLD,"endY = %d \t \n", endY);

	if(startY <= 1||endY >= n-1){
		PetscPrintf(PETSC_COMM_WORLD,"Wall collision! \n");
		return 1;
	}

	if(endX >= m-1){
		PetscPrintf(PETSC_COMM_WORLD,"Particle leaves the channel! \n");
		return 1;
	}

	double Fint = 0., Gint = 0., Mint = 0., Qint = 0.;
	double* Qmint = make1DDoubleArray(Ns);

	for(int i=startX; i<=endX; i++){
		double xV = (i-0.5)*h;
		for(int j=startY; j<=endY; j++){
			double yU = (j-0.5)*h;
			double f = -Ip_U[k][i][j]*(U[i][j]-Us[i][j])/dtau;
			double g = -Ip_V[k][i][j]*(V[i][j]-Vs[i][j])/dtau;
			double q = -Ip_S[k][i][j]*(T[i][j]-Ts[i][j])/dtau;
			double* qm = make1DDoubleArray(Ns);
			for(int s=0; s<Ns; s++){
				qm[s] = -Ip_S[k][i][j]*(C[s][i][j]-Cs[s][i][j])/dtau;
			}
			Fint += f*h*h; /* units : m/s^2 */
			Gint += g*h*h; /* units : m/s^2 */
			Mint +=((xV-xg[k])*g-(yU-yg[k])*f)*h*h;/* units: m^2/s^2 */
			Qint += q*h*h; /*units : K/s */
			for(int s=0; s<Ns; s++){
				Qmint[s] += qm[s]*h*h; /*units : mol/m.s */
			}
		}
	}

	dudt[k][2] = -Fint/(Sp[k]*(rho_r-1.));
	dvdt[k][2] = -Gint/(Sp[k]*(rho_r-1.));
	domegadt[k][2] = -Mint/(II[k]*Sp[k]*(rho_r-1.));
	for(int s=0; s<Ns; s++){
		dCdt[k][s][2] = -Qmint[s]/Sp[k]; /* units : mol/m3.s */ 
	}
	double Qr = dCdt[k][0][2]*Sp[k]*(-dH); /*units : J/m.s */ 
	dTdt[k][2] = -Qint/(Sp[k]*(rho_r*cr-1.)) + Qr/(Sp[k]*(rho_p*cp-rho_f*cf));

	/* Compute Hydrodynamic forces */ 

#ifdef MOVE
	F_drag[k] = rho_f*(-rho_r/(rho_r - 1.))*Fint; //[N/m]
	F_lift[k] = rho_f*(-rho_r/(rho_r - 1.))*Gint; //[N/m]
	Torque[k] = rho_f*(-rho_r/(rho_r - 1.))*Mint; //[N]
	Q_heat[k] = rho_f*cf*((-rho_r*cr/(rho_r*cr-1.))*Qint + Qr/(rho_p*cp-rho_f*cf)); //[W/m]
#endif

#ifndef MOVE
	F_drag[k] = -rho_f*Fint; // [N]
	F_lift[k] = -rho_f*Gint;
	Torque[k] = -rho_f*Mint; //[N.m]
	Q_heat[k] = -rho_f*cf*Qint; // [W]
	Phi_species[k] = -Qmint[0]; //[mol/s]
	// Here, Phi_species corresponds to the flux of A (reactant) at the surface of particle 0 //
#endif

	PetscPrintf(PETSC_COMM_WORLD,"Force along -x dir on particle %d = %1.6e [N/m]  \n", k+1, F_drag[k]);
	PetscPrintf(PETSC_COMM_WORLD,"Force along -y dir on particle %d = %1.6e [N/m]  \n", k+1, F_lift[k]);
	PetscPrintf(PETSC_COMM_WORLD,"Torque on particle %d = %1.6e [N]  \n", k+1, Torque[k]);
	PetscPrintf(PETSC_COMM_WORLD,"Heat flux on particle %d = %1.6e [W/m] \n", k+1, Q_heat[k]);
	return 0;
}

void get_ghosts(double** U, double** V, double** P, double** Told, double*** Cold, double CA0, double CB0)
{
	/*Ghost points to impose BC's */
	for (int i=0; i<m; i++){
		/* On U */
		/* Bottom Wall : slip : du/dn= 0 */
		U[i][0] = U[i][1];
		/* Top Wall : slip : du/dn = 0  */
		U[i][n-1] = U[i][n-2];

		/* On P, T, C */
		/* Walls : dpdn = 0 */
		P[i][0] = P[i][1];
		P[i][n-1] = P[i][n-2];

		/* Walls : adiabatic: dTdn = 0, no mass flux */
		Told[i][0] = Told[i][1];
		Told[i][n-1] = Told[i][n-2];

		for (int s=0; s<Ns; s++){
			Cold[s][i][0] = Cold[s][i][1];
			Cold[s][i][n-1] = Cold[s][i][n-2];
		}
	}
	for (int j = 0; j<n; j++){
		/* On V */
		/* Inflow : horizontal flow --> V = 0 */
		V[0][j] = -0.2*(V[3][j] - 5.*V[2][j] + 15.*V[1][j]);
		/* Natural Ouflow : dV/dx =0 */
		V[m-1][j] = V[m-2][j];

		/* On P */
		/* Outflow : P = 0 */
		P[m-1][j] = -P[m-2][j];

		/* On T and C */
		/* Inflow : T uniform  */ 
		Told[0][j] = -0.2*(Told[3][j]-5.*Told[2][j]+15.*Told[1][j]-16.*Tm0);

		/* Inflow : CA = CA0; CB = CB0 */
		Cold[0][0][j] = -0.2*(Cold[0][3][j]-5.*Cold[0][2][j]+15.*Cold[0][1][j]-16.*CA0);
		Cold[1][0][j] = -0.2*(Cold[1][3][j]-5.*Cold[1][2][j]+15.*Cold[1][1][j]-16.*CB0);

		/*Outflow : We cancel axial dispersion d2T/dx2 = 0; d2C/dx2 = 0; */
		Told[m-1][j] = (7.*Told[m-2][j]-5.*Told[m-3][j]+Told[m-4][j])/3.;
		for(int s=0; s<Ns; s++){
			Cold[s][m-1][j] = (7.*Cold[s][m-2][j]-5.*Cold[s][m-3][j]+Cold[s][m-4][j])/3.;
		}
	}
	PetscPrintf(PETSC_COMM_WORLD, "CA ghost = %f [mol/m^3]\n", Cold[0][0][50]);
}

void get_masks(double*** Ip_S, double*** Ip_U, double*** Ip_V, double** I_S, double** I_U, double** I_V,  double* xg, double* yg, double* rp, double* theta, double** coloring)
{
	for(int i=0; i<m; i++){
		double xU = i*h;
		double xV = (i-0.5)*h;
		for(int j=0; j<n; j++){
			double yU = (j-0.5)*h;
			double yV = j*h;
			I_S[i][j] = 0; /*Reset the masks */
			I_U[i][j] = 0;
			I_V[i][j] = 0;
			coloring[i][j] = 0;

			/*Go over all the particles */
			for(int k=0; k<Np; k++){
				Ip_S[k][i][j]=((xV-xg[k])*(xV-xg[k])+(yU-yg[k])*(yU-yg[k])<= rp[k]*rp[k]);
				Ip_U[k][i][j]=((xU-xg[k])*(xU-xg[k])+(yU-yg[k])*(yU-yg[k])<= rp[k]*rp[k]);
				Ip_V[k][i][j]=((xV-xg[k])*(xV-xg[k])+(yV-yg[k])*(yV-yg[k])<= rp[k]*rp[k]);

				double xloc = xV-xg[k];
				double yloc = yU-yg[k];
				double delta = 0;
				coloring[i][j] += Ip_S[k][i][j];
				if (xloc>0 && yloc>0){ //1st Quadrant
					delta = atan(yloc/xloc);
				}
				if (xloc>0 && yloc<0){ //4th Quadrant
					delta = atan(yloc/xloc) + 2.*M_PI;
				}
				if (xloc<0){
					delta = atan(yloc/xloc) + M_PI;
				}
				if (xloc == 0 && yloc > 0){
					delta = M_PI/2.;
				}
				if (xloc == 0 && yloc < 0){
					delta = 3.*M_PI/2.;
				}

				if((int)((delta-theta[k])/(M_PI/2.)) % 2 == 0 ){
					coloring[i][j] = -coloring[i][j];
				}
				I_S[i][j] += Ip_S[k][i][j];
				I_U[i][j] += Ip_U[k][i][j];
				I_V[i][j] += Ip_V[k][i][j];
			}

			if(I_S[i][j] > 1 || I_U[i][j] > 1 || I_V[i][j] > 1 ){
				PetscPrintf(PETSC_COMM_WORLD, "Erreur : collision de particules \n");
			}
		}
	}
}

void get_Cs(double*** Cs, double*** Ip_S, double** Cp)
{
	for(int i=0; i<m; i++){
		for(int j=0; j<n; j++){
			Cs[1][i][j] = 0.;
			for(int k=0; k<Np; k++){
				Cs[1][i][j] += Ip_S[k][i][j]*Cp[k][1];
			}
		}
	}
}
void get_Ts(double** Ts, double*** Ip_S, double* Tp){
	for(int i=0; i<m; i++){
		for(int j=0; j<n; j++){
			Ts[i][j] = 0.;
			for(int k = 0; k<Np; k++){
				Ts[i][j]+= Ip_S[k][i][j]*Tp[k];

			}
		}
	}
}

void get_Us_Vs(double** Us, double** Vs, double*** Ip_U, double*** Ip_V, double** Up, double** Vp, double** Wp, double* xg, double* yg){
	for(int i=0; i<m; i++){
		double xV = (i-0.5)*h;
		for(int j=0; j<n; j++){
			double yU = (j-0.5)*h;
			Us[i][j] = 0.;
			Vs[i][j] = 0.;
			for (int k = 0; k<Np; k++){
				Us[i][j]+= Ip_U[k][i][j]*(Up[k][3] - Wp[k][3]*(yU-yg[k]));
				Vs[i][j]+= Ip_V[k][i][j]*(Vp[k][3] + Wp[k][3]*(xV-xg[k]));
			}
		}
	}
}


void get_Ustar(double** U, double** Ustar, double** V, double** A, double** Aold, double** Us, double** P, double** I_U, double ramp)
{
	/* Ustar ADAMS-BASHFORTH 2 */
	for (int i=1; i<m-2; i++){
		for (int j=1; j<n-1; j++){
			/*only interior points*/
			/*Ustar */
			double Vij = .25*(V[i][j]+V[i+1][j]+V[i][j-1]+V[i+1][j-1]); /* average of the four neighbour points */
			A[i][j]= U[i][j]*(U[i+1][j]-U[i-1][j])/(2.*h) + Vij*(U[i][j+1]-U[i][j-1])/(2.*h);
			Ustar[i][j] = (U[i][j] + dt*(-1.5*A[i][j]+0.5*Aold[i][j] - (P[i+1][j]-P[i][j])/h
						+ nu*(U[i+1][j]+U[i-1][j]+U[i][j+1]+U[i][j-1]-4.*U[i][j])/(h*h)
						+ ramp*I_U[i][j]*Us[i][j]/dtau))/(1.+ramp*I_U[i][j]*dt/dtau);

			Aold[i][j]=A[i][j];
		}
	}
	/*Outflow condition */
	for (int j=1; j<n-1; j++){
		Ustar[m-2][j] = U[m-2][j] - dt*Um*(U[m-2][j]-U[m-3][j])/h;
	}
	/*Ustar[0][j] (inflow) is fixed once for all at the beginning */
}

void get_Ustar_EE(double** U, double** Ustar, double** V, double** Aold, double** Us, double** P, double** I_U, double ramp)
{
	/* Ustar EULER EXPLICIT */\
		for (int i=1; i<m-2; i++){
			for (int j=1; j<n-1; j++){
				/*only interior points*/
				/*Ustar */
				double Vij = .25*(V[i][j]+V[i+1][j]+V[i][j-1]+V[i+1][j-1]); /* average of the four neighbour points */
				Aold[i][j]= U[i][j]*(U[i+1][j]-U[i-1][j])/(2.*h) + Vij*(U[i][j+1]-U[i][j-1])/(2.*h);
				Ustar[i][j] = (U[i][j] + dt*(-Aold[i][j] - (P[i+1][j]-P[i][j])/h
							+ nu*(U[i+1][j]+U[i-1][j]+U[i][j+1]+U[i][j-1]-4.*U[i][j])/(h*h)
							+ ramp*I_U[i][j]*Us[i][j]/dtau))/(1.+ramp*I_U[i][j]*dt/dtau);
			}
		}
	/*Outflow condition */
	for (int j=1; j<n-1; j++){
		Ustar[m-2][j] = U[m-2][j] - dt*Um*(U[m-2][j]-U[m-3][j])/h;
	}
	/*Ustar[0][j] (inflow) is fixed once for all at the beginning --> but it does not matter with the MAC mesh !*/
}

void get_Vstar(double** V, double** Vstar, double** U, double** B, double** Bold, double** Vs, double** P, double** I_V, double ramp)
{
	/* Vstar  ADAMS-BASHFORTH 2 */
	for (int i=1; i<m-1; i++){
		for (int j=1; j<n-2; j++){
			double Uij = .25*(U[i][j]+U[i-1][j]+U[i][j+1]+U[i-1][j+1]); /* average of the four neighbour points */
			B[i][j] = Uij*(V[i+1][j]-V[i-1][j])/(2.*h) + V[i][j]*(V[i][j+1]-V[i][j-1])/(2.*h);
			Vstar[i][j] = (V[i][j] + dt*(-1.5*B[i][j]+0.5*Bold[i][j] - (P[i][j+1]-P[i][j])/h
						+ nu*(V[i+1][j]+V[i-1][j]+V[i][j+1]+V[i][j-1]-4.*V[i][j])/(h*h)
						+ ramp*I_V[i][j]*Vs[i][j]/dtau))/(1.+ramp*I_V[i][j]*dt/dtau);
			Bold[i][j]=B[i][j];
			/* the value of Vstar on the boundaries (j=0, j=n-2) does not matter with the MAC mesh */
		}
	}
}

void get_Vstar_EE(double** V, double** Vstar, double** U, double** Bold, double** Vs, double** P, double** I_V, double ramp)
{
	/* Vstar EULER EXPLICIT */
	for (int i=1; i<m-1; i++){
		for (int j=1; j<n-2; j++){
			double Uij = .25*(U[i][j]+U[i-1][j]+U[i][j+1]+U[i-1][j+1]); /* average of the four neighbour points */
			Bold[i][j] = Uij*(V[i+1][j]-V[i-1][j])/(2.*h) + V[i][j]*(V[i][j+1]-V[i][j-1])/(2.*h);
			Vstar[i][j] = (V[i][j] + dt*(-Bold[i][j] - (P[i][j+1]-P[i][j])/h
						+ nu*(V[i+1][j]+V[i-1][j]+V[i][j+1]+V[i][j-1]-4.*V[i][j])/(h*h)
						+ ramp*I_V[i][j]*Vs[i][j]/dtau))/(1.+ramp*I_V[i][j]*dt/dtau);
		}
	}
}

void poisson_solver(double** Ustar, double** Vstar, double **phi, int myrank, int nbproc)
{
	/* Solve the linear system Ax = b for a 2-D poisson equation on a structured grid */
	KSP sles;
	Mat A;
	Vec b, x;

	double div_u_star;
	int r, rowStart, rowEnd, i, j, ii, jj, its;
	int mytag = 12;
	int my_rowStart, my_rowEnd;
	double*  my_array;
	double* array = malloc(M*N*sizeof(double));
	MPI_Status status[3];

	/* Create the Laplacian matrix : A  */
	MatCreate( PETSC_COMM_WORLD, &A );
	MatSetSizes(A, PETSC_DECIDE, PETSC_DECIDE, M*N, M*N);
	//MatSetType(A, MATMPIAIJ);
	//MatMPIAIJSetPreallocation(A, 1, NULL, 4, NULL);
	MatSetUp(A);
	//MatSetFromOptions(A);
	//MatSetOption(A, MAT_NEW_NONZERO_LOCATION_ERR, PETSC_FALSE);
	MatGetOwnershipRange(A, &rowStart, &rowEnd);
	//MatAIJSetPreallocation(A, 5, NULL);

	//PetscPrintf(PETSC_COMM_WORLD, "starting row is %d \n", rowStart);
	//printf("starting row is %d \n", rowStart);
	for(r = rowStart; r<rowEnd; r++){
		ii = r/N; jj=r%N;
		if(ii>0){
			MatSetValue(A, r, r-N, -1., INSERT_VALUES);
		}
		if(jj>0){
			MatSetValue(A, r, r-1, -1., INSERT_VALUES);
		}
		MatSetValue(A, r, r, 4., INSERT_VALUES);
		if(jj<N-1){
			MatSetValue(A, r, r+1, -1., INSERT_VALUES);
		}
		if(ii<M-1){
			MatSetValue(A, r, r+N, -1., INSERT_VALUES);
		}
		if(ii == 0 || jj == 0 || jj == N-1)
		{
			MatSetValue(A, r, r, 3., INSERT_VALUES);
		}
		if(ii == 0 && (jj==0 || jj == N-1) )
		{
			MatSetValue(A, r, r, 2., INSERT_VALUES);
		}
		if(ii == M-1)
		{
			MatSetValue(A, r, r, 5., INSERT_VALUES);
		}
		if(ii == M-1 && (jj==0 || jj==N-1))
		{
			MatSetValue(A, r, r, 4., INSERT_VALUES);
		}
	}
	MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);

	/* Create the right-hand-side vector : b */
	VecCreate (PETSC_COMM_WORLD, &b );
	VecSetSizes(b, PETSC_DECIDE, M*N);
	VecSetFromOptions(b);
	VecGetOwnershipRange( b, &rowStart, &rowEnd );
	for(r = rowStart; r< rowEnd; r++){
		ii = r/N; jj=r%N;
		i = ii+1; j = jj+1;
		div_u_star = (Ustar[i][j]-Ustar[i-1][j]+Vstar[i][j]-Vstar[i][j-1])/h;

		VecSetValue(b, r, -(h*h/dt)*div_u_star, INSERT_VALUES);
	}
	VecAssemblyBegin(b);
	VecAssemblyEnd(b);

	/*Solve the linear system of equations */
	VecDuplicate( b, &x );
	KSPCreate(PETSC_COMM_WORLD, &sles );
	KSPSetOperators( sles, A, A );
	KSPSetFromOptions( sles );
	PetscPrintf(PETSC_COMM_WORLD,"Assembly is done \n");
	KSPSolve( sles, b, x );
	KSPGetIterationNumber( sles, &its );
	PetscPrintf( PETSC_COMM_WORLD, "Solution to Poisson eqn in %d iterations \n ", its);

	/*Transfer the solution to phi[i][j] */
	VecGetArray(x, &my_array);

	if ((M*N) % nbproc == 0){
		int length = ((M*N)/nbproc);
		MPI_Allgather(my_array, length, MPI_DOUBLE, array, length, MPI_DOUBLE, PETSC_COMM_WORLD);
		for (r = 0; r<M*N; r++){
			i = r/N; j = r%N;
			ii = i+1; jj = j+1;
			phi[ii][jj] = array[r];
		}
		for(ii=0; ii<m; ii++){
			/* cancel gradients : dp/dn=0 --> dphi/dn = 0*/
			phi[ii][0] = phi[ii][1];
			phi[ii][n-1] = phi[ii][n-2];
		}
		for(jj=0; jj<n; jj++){
			/*inflow : continuity of pressure gradient  */
			phi[0][jj] = phi[1][jj];
			/*outflow : zero pressure at outlet */
			phi[m-1][jj] = -phi[m-2][jj];
		}
	}
	else{
		if (myrank == 0){
			for (r=rowStart; r<rowEnd; r++) {
				array[r] = my_array[r];
			}
			for (int k = 1; k < nbproc; k++){
				MPI_Recv(&my_rowStart, 1, MPI_INT, k, mytag+1, PETSC_COMM_WORLD, &status[1] );
				MPI_Recv(&my_rowEnd, 1, MPI_INT, k, mytag+2, PETSC_COMM_WORLD, &status[2] );
				int length_k = my_rowEnd - my_rowStart;
				MPI_Recv(my_array, length_k, MPI_DOUBLE, k, mytag, PETSC_COMM_WORLD, &status[0] ) ;
				for (r=0; r<length_k; r++){
					int R = r+my_rowStart;
					array[R] = my_array[r];
				}
			}
		}
		else{
			MPI_Send(&rowStart, 1, MPI_INT, 0, mytag+1, PETSC_COMM_WORLD);
			MPI_Send(&rowEnd, 1, MPI_INT, 0, mytag+2, PETSC_COMM_WORLD);
			int length = rowEnd-rowStart;
			MPI_Send(my_array, length, MPI_DOUBLE, 0, mytag, PETSC_COMM_WORLD);
		}
		MPI_Bcast(array, M*N, MPI_DOUBLE, 0, PETSC_COMM_WORLD);

		for(r = 0; r< M*N; r++){
			i = r/N; j = r%N;
			ii = i+1; jj = j+1;
			phi[ii][jj] = array[r];
		}
		for(ii=0; ii<m; ii++){
			/* cancel gradients : dp/dn=0 --> dphi/dn = 0*/
			phi[ii][0] = phi[ii][1];
			phi[ii][n-1] = phi[ii][n-2];
		}
		for(jj=0; jj<n; jj++){
			/*inflow : continuity of pressure gradient  */
			phi[0][jj] = phi[1][jj];
			/*outflow : zero pressure at outlet */
			phi[m-1][jj] = -phi[m-2][jj];
		}
	}
	VecRestoreArray(x, &my_array);
	free(my_array);
	free(array);
	MatDestroy( &A );
	VecDestroy( &b ); VecDestroy( &x );
	KSPDestroy( &sles );

}

void old_poisson_solver(double** Ustar, double** Vstar, double **phi, double** R)
{
	double e = 0.;
	int SORiter = 0;
	do{
		for (int i=1; i<m-1; i++){
			for (int j=1; j<n-1; j++){
				double phistar = .25*(-(h/dt)*(Ustar[i][j]-Ustar[i-1][j]+Vstar[i][j]-Vstar[i][j-1])
						+ phi[i+1][j] + phi[i-1][j] + phi[i][j+1] + phi[i][j-1]);
				phi[i][j] = alpha*phistar + (1.-alpha)*phi[i][j];
			}
		}
		/*Adapt ghost points on Phi */
		for (int i=1; i<m-1; i++){
			phi[i][0] = phi[i][1];
			phi[i][n-1] = phi[i][n-2];
		}

		/*LEFT AND RIGHT */
		for (int j=1; j<n-1; j++){
			phi[0][j] = phi[1][j];
			phi[m-1][j] = -phi[m-2][j];
		}

		double Sum = 0.;

		for (int i=1; i<m-1; i++){
			for (int j=1; j<n-1; j++){
				R[i][j] = -(Ustar[i][j]-Ustar[i-1][j]+Vstar[i][j]-Vstar[i][j-1])/(h*dt)
					+(phi[i+1][j]+phi[i-1][j]+phi[i][j+1]+phi[i][j-1]-4.*phi[i][j])/(h*h);
				/* CACLULER LA SOMME */
				Sum += R[i][j]*R[i][j];

			}
		}

		e = (h*dt*d/Um)*pow(Sum/(L*d),0.5);
		SORiter ++;
	} while (e > SORtol && SORiter < SORitermax);
	PetscPrintf(PETSC_COMM_WORLD, "e = %1.3e after %d/%d iterations \n",e,SORiter,SORitermax);

}

void reset_phi(double** phi)
{
	for(int i=0; i<m; i++){
		for(int j=0; j<n; j++){
			phi[i][j] = 0.;
		}
	}
}


void update_flow(double** U, double** V, double** P, double** Ustar, double** Vstar, double** phi)
{
	for (int i=1; i<m-1; i++){
		int j;
		for (j=1; j<n-2; j++){
			U[i][j] = Ustar[i][j] - dt*(phi[i+1][j]-phi[i][j])/h;
			P[i][j] += phi[i][j];
			V[i][j] = Vstar[i][j] - dt*(phi[i][j+1]-phi[i][j])/h;
		}
		/* Last value of j is n-2 */
		U[i][j] = Ustar[i][j] - dt*(phi[i+1][j]-phi[i][j])/h;
		P[i][j] += phi[i][j];

		V[i][0] = (double) (4.*V[i][1]-V[i][2])/3.;/* to cancel gradient at the slip boundary : dv/dn = 0 */  /* Second-order accuracy */
		V[i][n-2] =(double) (4.*V[i][n-3]-V[i][n-4])/3.; /* to cancel gradient at the slip boundary : dv/dn = 0 */
	}
}

void update_temp_species(double** U, double** V,  double** T, double** Told, double** Ts, double** C_T, double** C_Told, double*** C, double*** Cold, double*** Cs, double*** C_C, double*** C_Cold, double** I_S, double ramp)
{
	/* Adams-Bashfort 2 */
	for (int i=1; i<m-1; i++){
		for (int j=1; j<n-1; j++){
			/*only interior points*/
			/* T */
			double Uij = .5*(U[i][j]+U[i-1][j]); /* average of the two neighbour points along x-direction */
			double Vij = .5*(V[i][j]+V[i-1][j]);
			C_T[i][j]= Uij*(Told[i+1][j]-Told[i-1][j])/(2.*h) + Vij*(Told[i][j+1]-Told[i][j-1])/(2.*h);
			T[i][j] = (Told[i][j] + dt*(-1.5*C_T[i][j] + 0.5*C_Told[i][j]
						+ alpha_f*(Told[i+1][j]+Told[i-1][j]+Told[i][j+1]+Told[i][j-1]-4.*Told[i][j])/(h*h)
						+ ramp*I_S[i][j]*Ts[i][j]/dtau))/(1.+ramp*I_S[i][j]*dt/dtau);
			for (int s=0; s<Ns; s++){
				C_C[s][i][j] = Uij*(Cold[s][i+1][j]-Cold[s][i-1][j])/(2.*h) + Vij*(Cold[s][i][j+1]-Cold[s][i][j-1])/(2.*h);
				C[s][i][j] = (Cold[s][i][j] + dt*(-1.5*C_C[s][i][j] + 0.5*C_Cold[s][i][j]
							+ Df[s]*(Cold[s][i+1][j]+Cold[s][i-1][j]+Cold[s][i][j+1]+Cold[s][i][j-1]-4.*Cold[s][i][j])/(h*h)
							+ ramp*I_S[i][j]*Cs[s][i][j]/dtau))/(1.+ramp*I_S[i][j]*dt/dtau);
			}
		}
	}
	for (int i=1; i<m-1; i++){
		for (int j=1; j<n-1; j++){
			C_Told[i][j] = C_T[i][j];
			Told[i][j] = T[i][j];
			for (int s=0; s<Ns; s++){
				C_Cold[s][i][j] = C_C[s][i][j];
				Cold[s][i][j] = C[s][i][j];
			}
		}
	}
}

void update_temp_species_EE(double** U, double** V, double** T, double** Told, double** Ts, double** C_Told, double*** C, double*** Cold, double*** Cs, double*** C_Cold, double** I_S, double ramp)
{
	/* Euler Explicit */
	for (int i=1; i<m-1; i++){
		for (int j=1; j<n-1; j++){
			/*only interior points*/
			/* T */
			double Uij = .5*(U[i][j]+U[i-1][j]); /* average of the two neighbour points along x-direction */
			double Vij = .5*(V[i][j]+V[i-1][j]);
			C_Told[i][j]= Uij*(Told[i+1][j]-Told[i-1][j])/(2.*h) + Vij*(Told[i][j+1]-Told[i][j-1])/(2.*h);
			T[i][j] = (Told[i][j] + dt*(-C_Told[i][j]
						+ alpha_f*(Told[i+1][j]+Told[i-1][j]+Told[i][j+1]+Told[i][j-1]-4.*Told[i][j])/(h*h)
						+ ramp*I_S[i][j]*Ts[i][j]/dtau))/(1.+ramp*I_S[i][j]*dt/dtau);
			for (int s=0; s<Ns; s++){
				C_Cold[s][i][j] = Uij*(Cold[s][i+1][j]-Cold[s][i-1][j])/(2.*h) + Vij*(Cold[s][i][j+1]-Cold[s][i][j-1])/(2.*h);
				C[s][i][j] = (Cold[s][i][j] + dt*(-C_Cold[s][i][j]
							+ Df[s]*(Cold[s][i+1][j]+Cold[s][i-1][j]+Cold[s][i][j+1]+Cold[s][i][j-1]-4.*Cold[s][i][j])/(h*h)
							+ ramp*I_S[i][j]*Cs[s][i][j]/dtau))/(1.+ramp*I_S[i][j]*dt/dtau);
			}
		}
	}
	for (int i=1; i<m-1; i++){
		for (int j=1; j<n-1; j++){
			Told[i][j] = T[i][j];
			for (int s=0; s<Ns; s++){
				Cold[s][i][j] = C[s][i][j];
			}
		}
	}
}

void update_Xp(double* xg, double* yg, double* theta, double** Up, double** Vp, double** Wp,  int k)
{
	xg[k] += dt*(23.*Up[k][2]-16.*Up[k][1]+5.*Up[k][0])/12.;\
		 yg[k] += dt*(23.*Vp[k][2]-16.*Vp[k][1]+5.*Vp[k][0])/12.;\
		 theta[k] += dt*(23.*Wp[k][2]-16.*Wp[k][1]+5.*Wp[k][0])/12.;\
		 PetscPrintf(PETSC_COMM_WORLD,"Position of the center of mass of particle %d: (x,y) = (%f,%f) \n", k+1, xg[k], yg[k]);
}

void update_Up(double** Up, double** Vp, double** Wp, double** F, double** G, double** M, int k)
{
	Up[k][3] = Up[k][2] + dt*(23.*F[k][2]-16.*F[k][1]+5.*F[k][0])/12.;\
		   Vp[k][3] = Vp[k][2] + dt*(23.*G[k][2]-16.*G[k][1]+5.*G[k][0])/12.;\
		   Wp[k][3] = Wp[k][2] + dt*(23.*M[k][2]-16.*M[k][1]+5.*M[k][0])/12.;\

		   Up[k][0] = Up[k][1]; Up[k][1] = Up[k][2]; Up[k][2] = Up[k][3];\
		   Vp[k][0] = Vp[k][1]; Vp[k][1] = Vp[k][2]; Vp[k][2] = Vp[k][3];\
		   Wp[k][0] = Wp[k][1]; Wp[k][1] = Wp[k][2]; Wp[k][2] = Wp[k][3];\
}


void update_Tp(double* Tp, double** Qp, int k)
{
	Tp[k] += dt*(23.*Qp[k][2]-16.*Qp[k][1]+5.*Qp[k][0])/12.;\
		 PetscPrintf(PETSC_COMM_WORLD,"Temperature of particle %d: Tp = %f[K] \n", k+1, Tp[k]);
}

void update_Cp(double** Cp, double*** Qmp, int k)
{
	Cp[k][0] = 0.;
	Cp[k][1] += dt*(23.*(Qmp[k][0][2]+Qmp[k][1][2])-16.*(Qmp[k][0][1]+Qmp[k][1][1])+5.*(Qmp[k][0][0]+Qmp[k][1][0]))/12.;
	//printf("Concentration of B in particle %d: Cp_B = %3.13e [mol/m3] \n", k+1, Cp[k][1]);
}

void writeFile(FILE* file, double **data, int iStart, int iEnd, int jStart, int jEnd)
{
	if (file != NULL){
		for (int i = iStart; i<iEnd; i++){
			int j;
			for (j = jStart; j<jEnd-1; j++) {
				fprintf(file, "%3.13e \t",data[i][j]);
			}
			fprintf(file, "%3.13e", data[i][j]);
			fprintf(file, "\n");
		}
	}
	else{
		printf(" An error occured in writeFile : invalid file pointer.\n");
	}
}

double* make1DDoubleArray(int arraySize) {
	double* theArray = (double*) calloc(arraySize, sizeof(double));
	return theArray;
}

double** make2DDoubleArray(int arraySizeX, int arraySizeY) {
	double** theArray;
	theArray = (double**) malloc(arraySizeX*sizeof(double*));
	for (int ix = 0; ix < arraySizeX; ix++){
		theArray[ix] =(double*) calloc(arraySizeY, sizeof(double));
	}
	return theArray;
}

double*** make3DDoubleArray(int numberOfparticles, int arraySizeX, int arraySizeY) {
	double*** theArray;
	theArray = (double***) malloc(numberOfparticles*sizeof(double**));
	for (int k = 0; k < numberOfparticles; k++) {
		theArray[k] = (double**) malloc(arraySizeX*sizeof(double*));
	}
	for (int k=0; k < numberOfparticles; k++) {
		for (int ix=0; ix < arraySizeX; ix++) {
			theArray[k][ix] = calloc(arraySizeY, sizeof(double));
		}
	}
	return theArray;
}

int** make2DIntArray(int arraySizeX, int arraySizeY) {
	int** theArray;
	theArray = (int**) malloc(arraySizeX*sizeof(int*));
	for (int ix = 0; ix < arraySizeX; ix++) {
		theArray[ix] =(int*) calloc(arraySizeY, sizeof(int));
	}
	return theArray;
}

int*** make3DIntArray(int numberOfparticles, int arraySizeX, int arraySizeY) {
	int*** theArray;
	theArray = (int***) malloc(numberOfparticles*sizeof(int**));
	for (int k = 0; k < numberOfparticles; k++) {
		theArray[k] = (int**) malloc(arraySizeX*sizeof(int*));
	}
	for (int k=0; k < numberOfparticles; k++) {
		for (int ix=0; ix < arraySizeX; ix++) {
			theArray[k][ix] = calloc(arraySizeY, sizeof(int));
		}
	}
	return theArray;
}
