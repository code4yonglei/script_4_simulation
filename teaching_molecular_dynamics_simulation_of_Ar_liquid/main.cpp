/*
	 Simple source code designed for teaching "Chemical Modelling"

	  Molecular Dynamics simulation of lennard-jones nano-particles

				Yong-Lei Wang @ Stockholm University

	Building instructions:
	   - c++ -O3 -o simu main.cpp
	   - ./simu < simu_cord256.input
	
													v1.0 - 2018-12-01
													v1.1 - 2019-01-12

													v2.0 - 2019-12-05
													v2.1 - 2022-12-14
													
	** Four files: head.h, main.cpp, simu_plot.py and simu_cord.input
	** Generate two files from simulation: simu_ener.txt and simu_traj.xyz
	** Output four PNG figures after running python script.

*/
#include "head.h"

double elapsed_time=0.0;
double dt=0.002;
double t;
int Natoms;
double *m=NULL;
double R=0.0, Rcut=2.5;
double R2cut=0.0, RIcut=0.0, RI2cut=0.0, RI6cut=0.0, RI12cut=0.0;
double **r=NULL, **v=NULL, **f=NULL;
double **r_old=NULL, **v_old=NULL, **f_old=NULL;
double lx=7.3660, ly=7.3660, lz=7.3660;
double Hlx=0.5*lx, Hly=0.5*ly, Hlz=0.5*lz;
double Vir=0.0, Vol=0.0;
double KB=1.0, T=0.0, P=0.0;
double tau=0.05; // thermostat time constant
double Tset=100; // units - K

// Variables used for energy and force calculations
double U=0.0, KE=0.0, TE=0.0;
double dx=0.0, dy=0.0, dz=0.0;
double r2=0.0, r6=0.0, Ir6=0.0;
double px=0.0, py=0.0, pz=0.0;
double fcut=0.0, fijx=0.0, fijy=0.0, fijz=0.0, F=0.0;

int main(int argc, char * argv[]){

	const char filename[]="simu_cord.input"; // Loading coordinate file
	Natoms = get_natoms(filename); // Determine the number of atoms
	init(Natoms, filename); // Initialize positions for all atoms

	R2cut   = Rcut*Rcut; // parameters for cutoffs
	RIcut   = 1/Rcut;
	RI2cut  = 1/R2cut;
	RI6cut  = RI2cut*RI2cut*RI2cut;
	RI12cut = RI6cut*RI6cut;

	Vol = lx*ly*lz;

	std::ofstream simFile,enerFile;
	simFile.open("simu_traj.xyz");
	enerFile.open("simu_ener.txt");
	std::cout << "-------------------------------------------------------------------\
----------------------------------------------------------------" << std::endl;
	std::cout << "Timestep" << std::setw(12) << "Time" << std::setw(12) 
		<< "T" << std::setw(12) << "P" << std::setw(12) 
		<< "PE" << std::setw(12) << "KE" << std::setw(12) 
		<< "TE" << std::setw(15)
		<< "Px" << std::setw(15) << "Py" << std::setw(15) << "Pz" << std::setw(15)
		<< std::endl;
	std::cout << "-------------------------------------------------------------------\
----------------------------------------------------------------" << std::endl;

	for(int k=0; k<=50000; k++) { //Loop over timesteps 
		elapsed_time = dt*double(k);
		if(k==0) // Cal pair-energy and forces
			calc_virial_force();
		calc_kenergy();
		TE=U+KE;
		calc_inst_temp_pr();
		if(k%100 == 0) {
			std::cout << std::setw(8) << k 
			<< std::setw(12) << elapsed_time
			<< std::setw(15) << T  << std::setw(15) << P
			<< std::setw(12) << U  << std::setw(12) << KE
			<< std::setw(12) << TE
			<< std::setw(15) << px << std::setw(15) << py
			<< std::setw(15) << pz
			<< std::setw(15) << std::endl;
			calc_momentum();
			write_xyz(simFile, k);
			dump_stats(enerFile, k);
		}
		vv_scheme_nose_hoover(); // integration using NH thermostat
	}

	simFile.close();
	enerFile.close();

	return 0;
}







// ------------------------------------------------------------------ //
// ------------------------------------------------------------------ //

void init(int Natoms, const char* filename){

	srand (time(NULL));
	// Memory Management
	m     = new double[Natoms];  
	r     = new double*[Natoms];
	v     = new double*[Natoms];
	f     = new double*[Natoms];
	r_old = new double*[Natoms];
	v_old = new double*[Natoms];
	f_old = new double*[Natoms];

	for(int i=0; i<Natoms; i++){
		r[i]     = new double[3]; // rx, ry, rz
		v[i]     = new double[3]; // vx, vy, vz
		f[i]     = new double[3]; // fx, fy, fz
		r_old[i] = new double[3];
		v_old[i] = new double[3];
		f_old[i] = new double[3];
	}

	read_xyz(Natoms, filename); // initialize the positions
	for(int i=0; i<Natoms; i++){ // intialize mass, velocities, forces
		m[i] = 1.0;
		v[i][0] = 2*(double(rand())/double(RAND_MAX)) - 1.0;
		v[i][1] = 2*(double(rand())/double(RAND_MAX)) - 1.0;
		v[i][2] = 2*(double(rand())/double(RAND_MAX)) - 1.0;
		f[i][0] = 0.0;
		f[i][1] = 0.0;
		f[i][2] = 0.0;
	}

	// Normalize COM velocities so that total momentum is zero
	double Vx=0.0, Vy=0.0, Vz=0.0; // 
	for(int i=0; i < Natoms; i++){
	Vx += v[i][0];
	Vy += v[i][1];
	Vz += v[i][2];
	}
	for(int i=0; i<Natoms; i++){
	v[i][0] = v[i][0] - (1.0/double(Natoms))*Vx;
	v[i][1] = v[i][1] - (1.0/double(Natoms))*Vy;
	v[i][2] = v[i][2] - (1.0/double(Natoms))*Vz;
	}

	calc_kenergy(); // Cal temperature and re-scaling velocity
	calc_inst_temp_pr();
	double alpha = std::sqrt(Tset/T);
	for(int i=0; i<Natoms; i++){
		v[i][0] = alpha*v[i][0];
		v[i][1] = alpha*v[i][1];
		v[i][2] = alpha*v[i][2];
	}
}

// ------------------------------------------------------------------ //
// ------------------------------------------------------------------ //

void read_xyz(const int Natoms, const char* filename){
	int count=0;
	std::cout << "Loading coordinates..." << std::endl;
	std::ifstream in_xyz;
	in_xyz.open(filename, std::ifstream::in);
	while(count < Natoms){
		in_xyz >> r[count][0] >> r[count][1] >> r[count][2];
		count++;
	}
	in_xyz.close();
}

// ------------------------------------------------------------------ //
// ------------------------------------------------------------------ //

int get_natoms(const char * filename){
	int num_atom = 0;
	double rx,ry,rz; // few dummy variables
	std::ifstream in_xyz;
	in_xyz.open(filename, std::ifstream::in);
	while( !in_xyz.eof() ){
		in_xyz >> rx >> ry >> rz;
		num_atom++;
	}
	in_xyz.close();
	return (num_atom-1);
}

// ------------------------------------------------------------------ //
// ------------------------------------------------------------------ //

void calc_kenergy(void){
	KE = 0.0;
	for(int i=0; i<Natoms; i++)
		KE += (v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2]);
	KE *= 0.5;
}

// ------------------------------------------------------------------ //
// ------------------------------------------------------------------ //

void calc_inst_temp_pr(void){
	T=0.0;
	P=0.0; 
	T = (2*KE)/(3.0*(Natoms-1)*KB);
	P = ((Natoms*T)/Vol + (Vir/(3*Vol)))*42.49; // MPa
	T = T*120.962;
}

// ------------------------------------------------------------------ //
// ------------------------------------------------------------------ //

void calc_virial_force(void){

	dx=0.0; dy=0.0; dz=0.0; // rezero parameters
	r2 = 0.0; r6=0.0; Ir6=0.0;
	fcut=0.0; F=0.0; Vir=0.0;

	fcut = RIcut*(-48*RI12cut + 24*RI6cut);  // dU/dr at Rcut

	for(int i=0; i<Natoms; i++){ // forces acting on specific atoms
		f[i][0] = 0.0;
		f[i][1] = 0.0;
		f[i][2] = 0.0;
	}

	for(int i=0; i<Natoms; i++){
		for(int j=i+1; j<Natoms; j++){
			fijx = 0.0;
			fijy = 0.0;
			fijz = 0.0;
			R = 0;

			dx = r[i][0] - r[j][0];
			dy = r[i][1] - r[j][1];
			dz = r[i][2] - r[j][2];
		
			if(dx > Hlx) dx=dx-lx; // Apply PBC: Nearest Image Convention
			if(dx <-Hlx) dx=dx+lx;
			if(dy > Hly) dy=dy-ly;
			if(dy <-Hly) dy=dy+ly;
			if(dz > Hlz) dz=dz-lz;
			if(dz <-Hlz) dz=dz+lz;
			r2 = (dx*dx) + (dy*dy) + (dz*dz);


			if(r2<R2cut){
				r6 = r2*r2*r2;
				Ir6 = 1.0/r6;
				R = std::sqrt(r2);

				F    = 24*(2*Ir6*Ir6 - Ir6);
				fijx = F*(dx/r2) + fcut*(dx/R);
				fijy = F*(dy/r2) + fcut*(dy/R);
				fijz = F*(dz/r2) + fcut*(dz/R);
				
				f[i][0] += fijx;
				f[i][1] += fijy;
				f[i][2] += fijz;
				f[j][0] -= fijx;
				f[j][1] -= fijy;
				f[j][2] -= fijz;
				
				Vir += (F + (fcut*R));
			} // r2< R2cut
		}
	}
}

// ------------------------------------------------------------------ //
// ------------------------------------------------------------------ //

void vv_scheme_nose_hoover(void){
	double tfact;
	double eta_t=0.0, eta_tdt=0.0;
  
	for(int i=0; i<Natoms; i++){
		r_old[i][0] = r[i][0];
		r_old[i][1] = r[i][1];
		r_old[i][2] = r[i][2];
		v_old[i][0] = v[i][0];
		v_old[i][1] = v[i][1];
		v_old[i][2] = v[i][2];
		f_old[i][0] = f[i][0];
		f_old[i][1] = f[i][1];
		f_old[i][2] = f[i][2];
	}
	eta_t = (1/(tau*tau)) * ((T/Tset) -1);
    
	for(int i=0; i<Natoms; i++){
		tfact = dt/(2.0*m[i]);

		v[i][0] = v_old[i][0] + (f_old[i][0] - (eta_t*v_old[i][0]))*tfact; // Step 1. v(t+0.5dt)
		v[i][1] = v_old[i][1] + (f_old[i][1] - (eta_t*v_old[i][1]))*tfact;
		v[i][2] = v_old[i][2] + (f_old[i][2] - (eta_t*v_old[i][2]))*tfact;
		v_old[i][0] = v[i][0];
		v_old[i][1] = v[i][1];
		v_old[i][2] = v[i][2];

		r[i][0] = r_old[i][0] + v_old[i][0]*dt; // Step 2. r(t+dt)
		r[i][1] = r_old[i][1] + v_old[i][1]*dt;
		r[i][2] = r_old[i][2] + v_old[i][2]*dt;

		if(r[i][0] <  0) r[i][0] = r[i][0] + lx; // PBC
		if(r[i][0] > lx) r[i][0] = r[i][0] - lx;
		if(r[i][1] <  0) r[i][1] = r[i][1] + ly;
		if(r[i][1] > ly) r[i][1] = r[i][1] - ly;
		if(r[i][2] <  0) r[i][2] = r[i][2] + lz;
		if(r[i][2] > lz) r[i][2] = r[i][2] - lz;
	}

	calc_pairenergy();
	calc_virial_force();
	eta_tdt = eta_t + (dt/(tau*tau)) * ((T/Tset) -1);
	for(int i=0; i<Natoms; i++){
		v[i][0] = (v_old[i][0] + f[i][0]*tfact)/(1 + eta_tdt*tfact); // Step 3. v(t+dt) 
		v[i][1] = (v_old[i][1] + f[i][1]*tfact)/(1 + eta_tdt*tfact);
		v[i][2] = (v_old[i][2] + f[i][2]*tfact)/(1 + eta_tdt*tfact);
	}
}

// ------------------------------------------------------------------ //
// ------------------------------------------------------------------ //

void calc_pairenergy(void){

	double URcut=0.0; U=0.0; // rezero params before calculating energy
	URcut=4.0*(RI12cut - RI6cut);
    
	for(int i=0; i<Natoms; i++){
		for(int j=i+1; j<Natoms; j++){
			R=0.0;
			dx=0.0; dy=0.0; dz=0.0;
			r2=0.0; r6=0.0; Ir6=0.0;

			dx = r[i][0] - r[j][0];
			dy = r[i][1] - r[j][1];
			dz = r[i][2] - r[j][2];
			if(dx > Hlx) dx=dx-lx; // PBC
			if(dx <-Hlx) dx=dx+lx;
			if(dy > Hly) dy=dy-ly;
			if(dy <-Hly) dy=dy+ly;
			if(dz > Hlz) dz=dz-lz;
			if(dz <-Hlz) dz=dz+lz;
			r2 = (dx*dx) + (dy*dy) + (dz*dz);
      
			if(r2<R2cut){
				r6 = r2*r2*r2;
				Ir6 = 1/r6;
				R=std::sqrt(r2);
				U += (4*((Ir6*Ir6)-Ir6) - URcut
					- (R-Rcut)*RIcut*(-48*RI12cut+24*RI6cut));
			} // r2<R2cut
		}
	}
}

// ------------------------------------------------------------------ //
// ------------------------------------------------------------------ //

void calc_momentum(void){
	px=0.0; py=0.0; pz=0.0;
	for(int i=0; i<Natoms; i++){
		px += m[i]*v[i][0];
		py += m[i]*v[i][1];
		pz += m[i]*v[i][2];
	}
}

// ------------------------------------------------------------------ //
// ------------------------------------------------------------------ //

void write_xyz(std::ofstream & simFile, int config){
	double Xcom=0.0, Ycom=0.0, Zcom=0.0;
	simFile << "    " << Natoms+1 << std::endl;
	simFile << " i =  " << config << std::endl;
	for(int i=0; i<Natoms; i++) {
		simFile << "H" 
				<< std::setw(15) << r[i][0] 
				<< std::setw(20) << r[i][1] 
				<< std::setw(20) << r[i][2]
		<< std::endl;
		Xcom += r[i][0];
		Ycom += r[i][1];
		Zcom += r[i][2];
	}
	simFile << "O" // output COM as one atom
			<< std::setw(15) << Xcom/Natoms 
			<< std::setw(20) << Ycom/Natoms 
			<< std::setw(20) << Zcom/Natoms
	<< std::endl;
}

// ------------------------------------------------------------------ //
// ------------------------------------------------------------------ //

void dump_stats(std::ofstream & enerFile, int config){
	enerFile << std::setw(8) << config
		<< std::setw(12) << elapsed_time
		<< std::setw(12) << T  << std::setw(12) << P
		<< std::setw(12) << U  << std::setw(12) << KE
		<< std::setw(12) << TE << std::setw(15) << px
		<< std::setw(15) << py << std::setw(15) << py << std::endl;
}

// ------------------------------------------------------------------ //
// ------------------------------------------------------------------ //
