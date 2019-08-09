#include<iostream>
#include<cmath>
#include"simplebounce.h"

// integral by trapezoidal rule
double integral(const double *integrand, const double dr, const int n){
	double result = 0.;
	for(int i=0; i<n-1; i++){
		result += (integrand[i] + integrand[i+1])*dr/2.;
	}
	return result;
}


scalarfield::scalarfield(const int nphi_, const int n_, const int rmax_, const int dim_) {
	phi = new double[n_*nphi_];
	n = n_;
	nphi = nphi_;
	rmax = rmax_;
	dim = dim_;
	dr = rmax_ / (n_-1.);
	drinv = 1./dr;
	rinv = new double[n_];
	rinvCalc();
}
scalarfield::~scalarfield(){
	delete[] phi;
	delete[] rinv;
}
double scalarfield::r(const int i) const {
	return dr*i;
}
double scalarfield::val(const int i, const int iphi) const {
	return phi[i*nphi + iphi];
}
void scalarfield::set(const int i, const int iphi, const double phi_){
	phi[i*nphi + iphi] = phi_;
}

// \nabla^2 \phi = d^2 phi / dr^2 + (d-1)/r * dphi/dr
double scalarfield::lap(const int i, const int iphi) const {
	if(i==0){
		return 2.*(phi[1*nphi + iphi]-phi[0*nphi + iphi])/dr/dr* dim;
	} else {
		return (phi[(i+1)*nphi + iphi] - 2.*phi[i*nphi + iphi] + phi[(i-1)*nphi + iphi])*drinv*drinv
				+ (phi[(i+1)*nphi + iphi] - phi[(i-1)*nphi + iphi])*0.5*drinv * (dim-1.)*rinv[i];
	}
}

void scalarfield::rinvCalc(){
	for(int i=1; i<n; i++){
		rinv[i] = 1./r(i);
	}
}

bounce::bounce() : scalarfield(1, 100, 1., 4) {
	n = 100;
	rmax = 1.;
	dim = 4;
	nphi = 1;
	RHS = new double[n*nphi];
	phiTV = new double[nphi];
	phiFV = new double[nphi];
	r_dminusoneth = new double[n];
	for(int i=0; i<n; i++){
		r_dminusoneth[i] = pow(r(i),dim-1);
	}

	// flags
	setModelDone = false;
	setVacuumDone = false;
	verbose = true;

	// parameters for numerical calculation
	safetyfactor = 0.9;
	maximumvariation = 0.01;
	xTV0 = 0.5;
	width0 = 0.05;
	derivMax = 1e-2;
	tend0 = 0.05;
	tend1 = 0.4;
	maxN = 1000;
}

bounce::~bounce(){
	delete[] RHS;
	delete[] phiTV;
	delete[] phiFV;
	delete[] r_dminusoneth;
}

void bounce::setRmax(const double rmax_){
	rmax = rmax_;
	dr = rmax / (n-1.);
	drinv = 1./dr;
	rinvCalc();
	for(int i=0; i<n; i++){
		r_dminusoneth[i] = pow(r(i),dim-1);
	}
	if(verbose){
		std::cerr << std::endl;
		std::cerr << "maximum of radius is set."<< std::endl;
		std::cerr << "\trmax : " << rmax << std::endl;
		std::cerr << std::endl;
	}
}

void bounce::setDimension(const int dim_){
	dim = dim_;
	for(int i=0; i<n; i++){
		r_dminusoneth[i] = pow(r(i),dim-1);
	}
	if(verbose){
		std::cerr << std::endl;
		std::cerr << "dimension is set."<< std::endl;
		std::cerr << "\tdim : " << dim << std::endl;
		std::cerr << std::endl;
	}
}

// set the number of grid. grid spacing dr is consistently changed.
void bounce::setN(const int n_){
	n = n_;
	dr = rmax / (n_-1.);
	drinv = 1./dr;
	delete[] rinv;
	rinv = new double[n_];
	rinvCalc();
	delete[] phi;
	delete[] RHS;
	delete[] r_dminusoneth;
	phi = new double[n_*nphi];
	RHS = new double[n_*nphi];
	r_dminusoneth = new double[n_];
	for(int i=0; i<n_; i++){
		r_dminusoneth[i] = pow(r(i),dim-1);
	}

	if(verbose){
		std::cerr << std::endl;
		std::cerr << "number of grids is set."<< std::endl;
		std::cerr << "\t(n,dr) : (" << n << ", " << dr << ")" << std::endl;
		std::cerr << std::endl;
	}
}

// set scalar potential and its derivatives to be used for the bounce calculation
void bounce::setModel(genericModel * const model_){
	model = model_;
	nphi = model_->nphi;
	delete[] phi;
	delete[] phiTV;
	delete[] phiFV;
	delete[] RHS;
	phi = new double[n*nphi];
	phiTV = new double[nphi];
	phiFV = new double[nphi];
	RHS = new double[n*nphi];

	if(verbose){
		std::cerr << std::endl;
		std::cerr << "model has been set"<< std::endl;
		std::cerr << "\tnumber of fields :\t" << nphi << std::endl;
	}

	setModelDone = true;
}

// kinetic energy : \int_0^\infty dr r^{d-1} \sum_i (-1/2) \phi_i \nabla^2\phi_i
double bounce::t() const {
	double integrand[n];
	for(int i=0; i<n; i++){
		integrand[i] = 0.;
		for(int iphi=0; iphi<nphi; iphi++){
			integrand[i] += r_dminusoneth[i] * -0.5 * phi[i*nphi+iphi] * lap(i,iphi);
		}
	}
	return integral(integrand, dr, n);
}

// potential energy : \int_0^\infty dr r^{d-1} V(\phi)
double bounce::v() const{
	double integrand[n];
	for(int i=0; i<n; i++){
		integrand[i] = r_dminusoneth[i] * (model->vpot(&phi[i*nphi]) - VFV);
	}
	return integral(integrand, dr, n);
}

// evolve the configuration by ds
double bounce::evolve(const double ds){

	double laplacian[n*nphi];
	for(int i=0; i<n-1; i++){
		for(int iphi=0; iphi<nphi; iphi++){
			laplacian[i*nphi + iphi] = lap(i,iphi);
		}
	}

	// integral1 : \int_0^\infty dr r^{d-1} \sum_i (\partial V / \partial\phi_i) \nabla^2\phi_i
	// integral2 : \int_0^\infty dr r^{d-1} \sum_i (\partial V / \partial\phi_i)^2
	double integrand1[n], integrand2[n];
	//for(int i=0; i<n; i++){
	for(int i=0; i<n-1; i++){
		integrand1[i] = 0.;
		integrand2[i] = 0.;
		model->calcDvdphi(&phi[i*nphi]);
		for(int iphi=0; iphi<nphi; iphi++){
			// r_dminusoneth[i] is equal to pow(r(i),dim-1)
			integrand1[i] += r_dminusoneth[i] * model->dvdphi[iphi] * laplacian[i*nphi+iphi];
			integrand2[i] += r_dminusoneth[i] * model->dvdphi[iphi] * model->dvdphi[iphi];
		}
	}
	integrand1[n-1] = 0.;
	integrand2[n-1] = 0.;

	// Eq. 9 of 1907.02417
	lambda = integral(integrand1,dr,n) / integral(integrand2,dr,n);

	// RHS of Eq. 8 of 1907.02417
	// phi at boundary is fixed to phiFV and will not be updated.
	for(int i=0; i<n-1; i++){
		model->calcDvdphi(&phi[i*nphi]);
		for(int iphi=0; iphi<nphi; iphi++){
			RHS[i*nphi + iphi] = laplacian[i*nphi+iphi] - lambda*model->dvdphi[iphi];
		}
	}

	// if RHS of EOM at the origin is too big, smaller step is taken.
	double sum = 0.;
	for(int iphi=0; iphi<nphi; iphi++){
		sum += RHS[0*nphi+iphi]*RHS[0*nphi+iphi];
	}
	double dstilde = maximumvariation * fieldExcursion() / sqrt(sum);
	if(ds < dstilde){
		dstilde = ds;
	}

	// flow by Eq. 8 of 1907.02417
	// phi at boundary is fixed to phiFV and will not be updated.
	for(int i=0; i<n-1; i++){
		for(int iphi=0; iphi<nphi; iphi++){
			phi[i*nphi + iphi] += dstilde*RHS[i*nphi + iphi];
		}
	}

	return lambda;
}

// RHS of Eq. 8 of 1907.02417
double bounce::residual(const int i, const int iphi) const{
	model->calcDvdphi(&phi[i*nphi]);
	return lap(i,iphi) - lambda*model->dvdphi[iphi];
}

// RHS of EOM for the bounce solution
double bounce::residualBounce(const int i, const int iphi) const{
	model->calcDvdphi(&phi[i*nphi]);
	return lap(i,iphi)/lambda - model->dvdphi[iphi];
}

// Euclidean action in d-dimensional space 
double bounce::action() const{
	double area = dim * pow(M_PI,dim/2.) / tgamma(dim/2.+1.);
	double rescaled_t_plus_v = pow(lambda, dim/2.-1.)*t() + pow(lambda, dim/2.)*v();
	return area * rescaled_t_plus_v;
}


// boucne solution from scale transformation
double bounce::rBounce(const int i) const{
	return sqrt(lambda)*dr*i;
}

// set the posiiton of true and false vacua
int bounce::setVacuum(const double *phiTV_, const double *phiFV_){
	if(!setModelDone){
		std::cerr << std::endl;
		std::cerr << "!!! model has not been set yet !!!"<< std::endl;
		std::cerr << std::endl;
		return -1;
	}

	if(model->vpot(phiTV_) > model->vpot(phiFV_) ){
		std::cerr << "!!! energy of true vacuum is larger than false vacuum !!!" << std::endl;
		return -1;
	}
	for(int iphi=0; iphi<nphi; iphi++){
		phiTV[iphi] = phiTV_[iphi];
	}
	for(int iphi=0; iphi<nphi; iphi++){
		phiFV[iphi] = phiFV_[iphi];
	}

	if(verbose){
		std::cerr << std::endl;
		std::cerr << "true and false vacua have been set." << std::endl;

		std::cerr << "\tfalse vacuum : ";
		std::cerr << "(";
		for(int iphi=0; iphi<nphi-1; iphi++){
			std::cerr << phiFV_[iphi] << ", ";
		}
		std::cerr << phiFV_[nphi-1]<< ")\t";
		std::cerr << "V = " << model->vpot(phiFV_) << std::endl;

		std::cerr << "\ta point with smaller V : ";
		std::cerr << "(";
		for(int iphi=0; iphi<nphi-1; iphi++){
			std::cerr << phiTV_[iphi] << ", ";
		}
		std::cerr << phiTV_[nphi-1]<< ")\t";
		std::cerr << "V = " << model->vpot(phiTV_) << std::endl;
	}

	VFV = model->vpot(phiFV_);
	setVacuumDone = true;

	return 0;
};

// set the initial configuration
void bounce::setInitial(const double frac, const double width){
	//for(int i=0; i<n; i++){
	for(int i=0; i<n-1; i++){
		for(int iphi=0; iphi<nphi; iphi++){
			phi[i*nphi+iphi] = phiTV[iphi] + (phiFV[iphi]-phiTV[iphi])*(1.+tanh( (i-n*frac)/(n*width) ))/2.;
		}
	}
	for(int iphi=0; iphi<nphi; iphi++){
		phi[(n-1)*nphi+iphi] = phiFV[iphi];
	}
}

// field excursion from the origin to the infinity
double bounce::fieldExcursion() const{
	double normsquared = 0.;
	for(int iphi=0; iphi<nphi; iphi++){
		normsquared += pow(phi[(n-1)*nphi+iphi] - phi[0*nphi+iphi], 2);
	}
	return sqrt(normsquared);
}

// derivative of scalar field at boundary
double bounce::derivativeAtBoundary() const{
	double normsquared = 0.;
	for(int iphi=0; iphi<nphi; iphi++){
		normsquared += pow(phi[(n-1)*nphi+iphi] - phi[(n-2)*nphi+iphi], 2);
	}
	return sqrt(normsquared)/dr;
}

// evolve the configuration
double bounce::evolveUntil(const double tend){

	// 1 + d + sqrt(1 + d) is maximum of absolute value of eigenvalue of {{-2d, 2d},{ (1-d)/2 + 1, -2}},
	// which is discreitzed Laplacian for n = 2. This value is 6 for d=3, and 7.23607 for d=4.
	// The numerical value of maximum of absolute value of eigenvalue of discretized Laplacian for large n is 6 for d=3, and 7.21417 for d=4
	double t = 0.;
	double dt = 2./(1. + dim + sqrt(1.+dim)) * pow(dr,2) * safetyfactor;

	if(verbose){
		std::cerr << std::endl;
		std::cerr << "evolve until t = " << tend << ", (dt = " << dt << ")" << std::endl;
		std::cerr << std::endl;
	}

	while(t<tend){
		evolve(dt);
		t += dt;
	}
	return derivativeAtBoundary()/fieldExcursion();
}

// main routine to get the bounce solution
int bounce::solve(){
	if(!setModelDone){
		std::cerr << std::endl;
		std::cerr << "!!! model has not been set yet !!!"<< std::endl;
		std::cerr << std::endl;
		return -1;
	}
	if(!setVacuumDone){
		std::cerr << std::endl;
		std::cerr << "!!! correct vacua have not been set yet !!!"<< std::endl;
		std::cerr << std::endl;
		return -1;
	}


	// make the bubble wall thin to get negative potential energy 
	if(verbose){
		std::cerr << "========================================" << std::endl;
		std::cerr << std::endl;
		std::cerr << "probing a thickness to get negative V[phi] ..." << std::endl;
		std::cerr << std::endl;
	}
	double xTV = xTV0;
	double width = width0;
	while(true){
		setInitial(xTV, width);
		if(verbose){
			std::cerr << "\t" << "xTrueVacuum:\t" << xTV << std::endl;
			std::cerr << "\t" << "xWidth:\t" << width << std::endl;
			std::cerr << "\t" << "V[phi] :\t" << v() << std::endl;
			std::cerr << "\t" << "n :\t" << n << std::endl;
			std::cerr << "\t" << std::endl;
		}

		// OK if V is negative
		if(v() < 0.) {
			break;
		}

		// if V is positive, make the wall thin.
		width = width * 0.5;
		if(width*n < 1.) {
			if(verbose){
				std::cerr << std::endl;
				std::cerr << "the current mesh is too sparse. increase the number of points." << std::endl;
				std::cerr << std::endl;
			}
			setN(2*n);
		}

		if(n>maxN){
			std::cerr << "!!! n became too large !!!" << std::endl;
			return -1;
		}
	}

	// make the size of the bubble smaller enough than the size of the sphere
	if(verbose){
		std::cerr << std::endl;
		std::cerr << "probing the size of the bounce configuration ..." << std::endl;
		std::cerr << std::endl;
	}
	while(true){

		double deriv = evolveUntil(tend0);
		if(verbose){
			std::cerr << "\t" << "deriv :\t" << deriv << std::endl;
			std::cerr << "\t" << "field excursion :\t" << fieldExcursion() << std::endl;
			std::cerr << "\t" << "derivative at boundary:\t" << derivativeAtBoundary() << std::endl;
			std::cerr << "\t" << std::endl;
		}

		// dphi/dr at the boundary is small enough
		if( deriv  < derivMax) {
			break;
		// dphi/dr at the boundary is NOT small enough
		} else {

			// take smaller bounce configuration
			if(verbose){
				std::cerr << std::endl;
				std::cerr << "the size of the bounce is too large. initial condition is scale transformed." << std::endl;
				std::cerr << std::endl;
			}
			xTV = xTV * 0.5;
			width = width * 0.5;
			if(width*n < 1.) {
				if(verbose){
					std::cerr << std::endl;
					std::cerr << "the current mesh is too sparse. increase the number of points." << std::endl;
					std::cerr << std::endl;
				}
				setN(2*n);
			}

			// retry by using new initial condition
			setInitial(xTV, width);
			if(verbose){
				std::cerr << "\t" << "xTrueVacuum:\t" << xTV << std::endl;
				std::cerr << "\t" << "xWidth:\t" << width << std::endl;
				std::cerr << "\t" << "V[phi] :\t" << v() << std::endl;
				std::cerr << "\t" << "n :\t" << n << std::endl;
			}
		}

		if(n>maxN){
			std::cerr << "!!! n became too large !!!" << std::endl;
			return -1;
		}
	}

	if(verbose){
		std::cerr << std::endl;
		std::cerr << "minimizing the kinetic energy ..." << std::endl;
		std::cerr << std::endl;
	}

	evolveUntil(tend1);

	if(verbose){
		std::cerr << std::endl;
		std::cerr << "done." << std::endl;
		std::cerr << std::endl;
	}

	return 0;
}



double bounce::getlambda() const{
	return lambda;
}

// print the result
int bounce::printBounce() const{
	if(!setModelDone){
		std::cerr << std::endl;
		std::cerr << "!!! model has not been set yet !!!"<< std::endl;
		std::cerr << std::endl;
		return -1;
	}
	if(!setVacuumDone){
		std::cerr << std::endl;
		std::cerr << "!!! correct vacua have not been set yet !!!"<< std::endl;
		std::cerr << std::endl;
		return -1;
	}

	std::cout << "# ";
	std::cout << "r\t";
	for(int iphi=0; iphi<nphi; iphi++){
		std::cout << "phi["<< iphi <<"]\t";
	}
	for(int iphi=0; iphi<nphi; iphi++){
		std::cout << "RHS of EOM["<< iphi <<"]\t";
	}
	std::cout << std::endl;

	for(int i=0; i<n; i++){
		std::cout << rBounce(i) << "\t";
		for(int iphi=0; iphi<nphi; iphi++){
			std::cout << val(i,iphi) << "\t";
		}
		for(int iphi=0; iphi<nphi; iphi++){
			if(i!=(n-1)){
				std::cout << residualBounce(i,iphi) << "\t";
			} else {
				std::cout << 0. << "\t";
			}
		}
		std::cout << std::endl;
	}

	return 0;
}

