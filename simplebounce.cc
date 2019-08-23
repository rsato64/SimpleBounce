#include<iostream>
#include<cmath>
#include"simplebounce.h"

namespace simplebounce{

// integral by trapezoidal rule
double integral(const double *integrand, const double dr, const int n){
	double result = 0.;
	for(int i=0; i<n-1; i++){
		result += (integrand[i] + integrand[i+1])*dr*0.5;
	}
	return result;
}



// class Scalarfield

Scalarfield::Scalarfield(const int nphi__, const int n__, const int rmax__, const int dim__) {
	n_ = n__;
	nphi_ = nphi__;
	rmax_ = rmax__;
	dim_ = dim__;
	dr_ = rmax__ / (n__-1.);
	drinv_ = 1./dr_;
	phi_ = new double[n__*nphi__];
	rinv_ = new double[n__];
	r_dminusoneth_ = new double[n__];
	updateInfo();
}

Scalarfield::~Scalarfield(){
	delete[] phi_;
	delete[] rinv_;
	delete[] r_dminusoneth_;
}

// phi_iphi at r_i
double Scalarfield::phi(const int i, const int iphi) const {
	return phi_[i*nphi_+iphi];
}
void Scalarfield::setPhi(const int i, const int iphi, const double phi__) {
	phi_[i*nphi_+iphi] = phi__;
}
void Scalarfield::addToPhi(const int i, const int iphi, const double phi__) {
	phi_[i*nphi_+iphi] += phi__;
}
double* Scalarfield::phivec(const int i) const {
	return &phi_[i*nphi_+0];
}

// radius : r_i
double Scalarfield::r(const int i) const {
	return dr_*i;
}

// Laplacian in radial coordinate : \nabla^2 \phi = d^2 phi / dr^2 + (d-1)/r * dphi/dr
double Scalarfield::lap(const int i, const int iphi) const {
	if(i==0){
		return 2.*(phi(1,iphi)-phi(0,iphi))*drinv_*drinv_* dim_;
	} else {
		return (phi(i+1,iphi) - 2.*phi(i,iphi) + phi(i-1,iphi))*drinv_*drinv_
				+ (phi(i+1,iphi) - phi(i-1,iphi))*0.5*drinv_ * (dim_-1.)*rinv_[i];
	}
}

// calculate and save 1/r(i) and pow(r(i), dim-1)
void Scalarfield::updateInfo(){
	for(int i=1; i<n(); i++){
		rinv_[i] = 1./r(i);
	}
	for(int i=0; i<n(); i++){
		r_dminusoneth_[i] = pow(r(i),dim_-1);
	}
}

// set the radius at the boundary
void Scalarfield::setRmax(const double rmax__){
	rmax_ = rmax__;
	dr_ = rmax_ / (n_-1.);
	drinv_ = 1./dr_;
	updateInfo();
}

// set the dimension of the Euclidean space
void Scalarfield::setDimension(const int dim__){
	dim_ = dim__;
	updateInfo();
}

// set the number of grid. grid spacing dr is consistently changed.
void Scalarfield::setN(const int n__){
	n_ = n__;
	dr_ = rmax_ / (n__-1.);
	drinv_ = 1./dr_;
	delete[] phi_;
	delete[] rinv_;
	delete[] r_dminusoneth_;
	phi_ = new double[n__*nphi_];
	rinv_ = new double[n__];
	r_dminusoneth_ = new double[n__];
	updateInfo();
}

void Scalarfield::setNphi(const int nphi__){
	nphi_ = nphi__;
	delete[] phi_;
	phi_ = new double[n_*nphi_];
}

int Scalarfield::n() const{
	return n_;
}
int Scalarfield::nphi() const{
	return nphi_;
}
int Scalarfield::dim() const{
	return dim_;
}
double Scalarfield::rmax() const{
	return rmax_;
}
double Scalarfield::dr() const{
	return dr_;
}
double Scalarfield::r_dminusoneth(const int i) const{
	return r_dminusoneth_[i];
}


// class GenericModel
GenericModel::GenericModel(){
	dvdphi = new double[1];
}
GenericModel::~GenericModel(){
	delete[] dvdphi;
}
void GenericModel::setNphi(const int nphi_){
	nphi = nphi_;
	delete[] dvdphi;
	dvdphi = new double[nphi_];
}

// class BounceCalculator

BounceCalculator::BounceCalculator() : Scalarfield(1, 100, 1., 4) {
	phiTV = new double[nphi()];
	phiFV = new double[nphi()];

	// flags
	setModelDone = false;
	setVacuumDone = false;
	verbose = false;

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

BounceCalculator::~BounceCalculator(){
	delete[] phiTV;
	delete[] phiFV;
}

// set scalar potential and its derivatives to be used for the bounce calculation
void BounceCalculator::setModel(GenericModel * const model_){
	model = model_;
	setNphi(model_->nphi);
	delete[] phiTV;
	delete[] phiFV;
	phiTV = new double[nphi()];
	phiFV = new double[nphi()];

	if(verbose){
		std::cerr << "model has been set"<< std::endl;
		std::cerr << "\tnumber of fields :\t" << nphi() << std::endl;
	}

	setModelDone = true;
}

// kinetic energy : \int_0^\infty dr r^{d-1} \sum_i (-1/2) \phi_i \nabla^2\phi_i
double BounceCalculator::t() const {
	double integrand[n()];
	for(int i=0; i<n(); i++){
		integrand[i] = 0.;
		for(int iphi=0; iphi<nphi(); iphi++){
			integrand[i] += r_dminusoneth(i) * -0.5 * phi(i,iphi) * lap(i,iphi);
		}
	}
	return integral(integrand, dr(), n());
}

// potential energy : \int_0^\infty dr r^{d-1} V(\phi)
double BounceCalculator::v() const{
	double integrand[n()];
	for(int i=0; i<n(); i++){
		integrand[i] = r_dminusoneth(i) * (model->vpot(phivec(i)) - VFV);
	}
	return integral(integrand, dr(), n());
}

// evolve the configuration by dtau
double BounceCalculator::evolve(const double dtau){

	double laplacian[n()][nphi()];

	for(int i=0; i<n()-1; i++){
		for(int iphi=0; iphi<nphi(); iphi++){
			laplacian[i][iphi] = lap(i,iphi);
		}
	}

	// integral1 : \int_0^\infty dr r^{d-1} \sum_i (\partial V / \partial\phi_i) \nabla^2\phi_i
	// integral2 : \int_0^\infty dr r^{d-1} \sum_i (\partial V / \partial\phi_i)^2
	double integrand1[n()], integrand2[n()];
	for(int i=0; i<n()-1; i++){
		integrand1[i] = 0.;
		integrand2[i] = 0.;
		model->calcDvdphi(phivec(i));
		for(int iphi=0; iphi<nphi(); iphi++){
			integrand1[i] += r_dminusoneth(i) * model->dvdphi[iphi] * laplacian[i][iphi];
			integrand2[i] += r_dminusoneth(i) * model->dvdphi[iphi] * model->dvdphi[iphi];
		}
	}
	integrand1[n()-1] = 0.;
	integrand2[n()-1] = 0.;

	// Eq. 9 of 1907.02417
	lambda = integral(integrand1,dr(),n()) / integral(integrand2,dr(),n());

	double RHS[n()][nphi()];
	// RHS of Eq. 8 of 1907.02417
	// phi at boundary is fixed to phiFV and will not be updated.
	for(int i=0; i<n()-1; i++){
		model->calcDvdphi(phivec(i));
		for(int iphi=0; iphi<nphi(); iphi++){
			RHS[i][iphi] = laplacian[i][iphi] - lambda*model->dvdphi[iphi];
		}
	}

	// if RHS of EOM at the origin is too big, smaller step is taken.
	double sum = 0.;
	for(int iphi=0; iphi<nphi(); iphi++){
		sum += RHS[0][iphi]*RHS[0][iphi];
	}
	double dtautilde = maximumvariation * fieldExcursion() / sqrt(sum);
	if(dtau < dtautilde){
		dtautilde = dtau;
	}

	// flow by Eq. 8 of 1907.02417
	// phi at boundary is fixed to phiFV and will not be updated.
	for(int i=0; i<n()-1; i++){
		for(int iphi=0; iphi<nphi(); iphi++){
			addToPhi(i, iphi, dtautilde*RHS[i][iphi]);
		}
	}

	return lambda;
}

// RHS of Eq. 8 of 1907.02417
double BounceCalculator::residual(const int i, const int iphi) const{
	model->calcDvdphi(phivec(i));
	return lap(i,iphi) - lambda*model->dvdphi[iphi];
}

// RHS of EOM for the bounce solution
double BounceCalculator::residualBounce(const int i, const int iphi) const{
	model->calcDvdphi(phivec(i));
	return lap(i,iphi)/lambda - model->dvdphi[iphi];
}

// Euclidean action in d-dimensional space 
double BounceCalculator::action() const{
	double area = dim() * pow(M_PI,dim()/2.) / tgamma(dim()/2.+1.);
	double rescaled_t_plus_v = pow(lambda, dim()/2.-1.)*t() + pow(lambda, dim()/2.)*v();
	return area * rescaled_t_plus_v;
}

// boucne solution from scale transformation
double BounceCalculator::rBounce(const int i) const{
	return sqrt(lambda)*dr()*i;
}

// set the posiiton of true and false vacua
int BounceCalculator::setVacuum(const double *phiTV_, const double *phiFV_){
	if(!setModelDone){
		std::cerr << "!!! model has not been set yet !!!"<< std::endl;
		return -1;
	}

	if(model->vpot(phiTV_) > model->vpot(phiFV_) ){
		std::cerr << "!!! energy of true vacuum is larger than false vacuum !!!" << std::endl;
		return -1;
	}
	for(int iphi=0; iphi<nphi(); iphi++){
		phiTV[iphi] = phiTV_[iphi];
	}
	for(int iphi=0; iphi<nphi(); iphi++){
		phiFV[iphi] = phiFV_[iphi];
	}

	if(verbose){
		std::cerr << "true and false vacua have been set." << std::endl;

		std::cerr << "\tfalse vacuum : ";
		std::cerr << "(";
		for(int iphi=0; iphi<nphi()-1; iphi++){
			std::cerr << phiFV_[iphi] << ", ";
		}
		std::cerr << phiFV_[nphi()-1]<< ")\t";
		std::cerr << "V = " << model->vpot(phiFV_) << std::endl;

		std::cerr << "\ta point with smaller V : ";
		std::cerr << "(";
		for(int iphi=0; iphi<nphi()-1; iphi++){
			std::cerr << phiTV_[iphi] << ", ";
		}
		std::cerr << phiTV_[nphi()-1]<< ")\t";
		std::cerr << "V = " << model->vpot(phiTV_) << std::endl;
	}

	VFV = model->vpot(phiFV_);
	setVacuumDone = true;

	return 0;
};

// set the initial configuration
void BounceCalculator::setInitial(const double frac, const double width){
	for(int i=0; i<n()-1; i++){
		for(int iphi=0; iphi<nphi(); iphi++){
			setPhi(i, iphi, phiTV[iphi] + (phiFV[iphi]-phiTV[iphi])*(1.+tanh( (i-n()*frac)/(n()*width) ))/2. );
		}
	}
	for(int iphi=0; iphi<nphi(); iphi++){
		setPhi(n()-1, iphi, phiFV[iphi]);
	}
}

// field excursion from the origin to the infinity
double BounceCalculator::fieldExcursion() const{
	double normsquared = 0.;
	for(int iphi=0; iphi<nphi(); iphi++){
		normsquared += pow(phi(n()-1,iphi) - phi(0,iphi), 2);
	}
	return sqrt(normsquared);
}

// derivative of scalar field at boundary
double BounceCalculator::derivativeAtBoundary() const{
	double normsquared = 0.;
	for(int iphi=0; iphi<nphi(); iphi++){
		normsquared += pow(phi(n()-1,iphi) - phi(n()-2,iphi), 2);
	}
	return sqrt(normsquared)/dr();
}

// evolve the configuration
double BounceCalculator::evolveUntil(const double tauend){

	// 1 + d + sqrt(1 + d) is maximum of absolute value of eigenvalue of {{-2d, 2d},{ (1-d)/2 + 1, -2}},
	// which is discreitzed Laplacian for n = 2. This value is 6 for d=3, and 7.23607 for d=4.
	// The numerical value of maximum of absolute value of eigenvalue of discretized Laplacian for large n is 6 for d=3, and 7.21417 for d=4
	double tau = 0.;
	double dtau = 2./(1. + dim() + sqrt(1.+dim())) * pow(dr(),2) * safetyfactor;

	if(verbose){
		std::cerr << "evolve until tau = " << tauend << ", (dtau = " << dtau << ")" << std::endl;
	}

	while(tau<tauend){
		evolve(dtau);
		tau += dtau;
	}
	return derivativeAtBoundary()/fieldExcursion();
}

// main routine to get the bounce solution
int BounceCalculator::solve(){
	if(!setModelDone){
		std::cerr << "!!! model has not been set yet !!!"<< std::endl;
		return -1;
	}
	if(!setVacuumDone){
		std::cerr << "!!! correct vacua have not been set yet !!!"<< std::endl;
		return -1;
	}


	// make the bubble wall thin to get negative potential energy 
	if(verbose){
		std::cerr << "probing a thickness to get negative V[phi] ..." << std::endl;
	}
	double xTV = xTV0;
	double width = width0;
	while(true){
		setInitial(xTV, width);
		if(verbose){
			std::cerr << "\t" << "xTrueVacuum:\t" << xTV << std::endl;
			std::cerr << "\t" << "xWidth:\t" << width << std::endl;
			std::cerr << "\t" << "V[phi] :\t" << v() << std::endl;
			std::cerr << "\t" << "n :\t" << n() << std::endl;
		}

		// OK if V is negative
		if(v() < 0.) {
			break;
		}

		// if V is positive, make the wall thin.
		width = width * 0.5;
		if(width*n() < 1.) {
			if(verbose){
				std::cerr << "the current mesh is too sparse. increase the number of points." << std::endl;
			}
			setN(2*n());
		}

		if(n()>maxN){
			std::cerr << "!!! n became too large !!!" << std::endl;
			return -1;
		}
	}

	// make the size of the bubble smaller enough than the size of the sphere
	if(verbose){
		std::cerr << "probing the size of the bounce configuration ..." << std::endl;
	}
	while(true){

		double deriv = evolveUntil(tend0);
		if(verbose){
			std::cerr << "\t" << "deriv :\t" << deriv << std::endl;
			std::cerr << "\t" << "field excursion :\t" << fieldExcursion() << std::endl;
			std::cerr << "\t" << "derivative at boundary:\t" << derivativeAtBoundary() << std::endl;
		}

		// dphi/dr at the boundary is small enough
		if( deriv  < derivMax) {
			break;
		// dphi/dr at the boundary is NOT small enough
		} else {

			// take smaller bounce configuration
			if(verbose){
				std::cerr << "the size of the bounce is too large. initial condition is scale transformed." << std::endl;
			}
			xTV = xTV * 0.5;
			width = width * 0.5;
			if(width*n() < 1.) {
				if(verbose){
					std::cerr << "the current mesh is too sparse. increase the number of points." << std::endl;
				}
				setN(2*n());
			}

			// retry by using new initial condition
			setInitial(xTV, width);
			if(verbose){
				std::cerr << "\t" << "xTrueVacuum:\t" << xTV << std::endl;
				std::cerr << "\t" << "xWidth:\t" << width << std::endl;
				std::cerr << "\t" << "V[phi] :\t" << v() << std::endl;
				std::cerr << "\t" << "n :\t" << n() << std::endl;
			}
		}

		if(n()>maxN){
			std::cerr << "!!! n became too large !!!" << std::endl;
			return -1;
		}
	}

	if(verbose){
		std::cerr << "minimizing the kinetic energy ..." << std::endl;
	}

	evolveUntil(tend1);

	if(verbose){
		std::cerr << "done." << std::endl;
	}

	return 0;
}

double BounceCalculator::getlambda() const{
	return lambda;
}

// print the result
int BounceCalculator::printBounce() const{
	if(!setModelDone){
		std::cerr << "!!! model has not been set yet !!!"<< std::endl;
		return -1;
	}
	if(!setVacuumDone){
		std::cerr << "!!! correct vacua have not been set yet !!!"<< std::endl;
		return -1;
	}

	std::cout << "# ";
	std::cout << "r\t";
	for(int iphi=0; iphi<nphi(); iphi++){
		std::cout << "phi["<< iphi <<"]\t";
	}
	for(int iphi=0; iphi<nphi(); iphi++){
		std::cout << "RHS of EOM["<< iphi <<"]\t";
	}
	std::cout << std::endl;

	for(int i=0; i<n(); i++){
		std::cout << rBounce(i) << "\t";
		for(int iphi=0; iphi<nphi(); iphi++){
			std::cout << phi(i,iphi) << "\t";
		}
		for(int iphi=0; iphi<nphi(); iphi++){
			if(i!=(n()-1)){
				std::cout << residualBounce(i,iphi) << "\t";
			} else {
				std::cout << 0. << "\t";
			}
		}
		std::cout << std::endl;
	}

	return 0;
}

void BounceCalculator::setSafetyfactor(double x){
	safetyfactor = x;
}
void BounceCalculator::setMaximumvariation(double x){
	maximumvariation = x;
}
void BounceCalculator::setXTV0(double x){
	xTV0 = x;
}
void BounceCalculator::setWidth0(double x){
	width0 = x;
}
void BounceCalculator::setDerivMax(double x){
	derivMax = x;
}
void BounceCalculator::setTend0(double x){
	tend0 = x;
}
void BounceCalculator::setTend1(double x){
	tend1 = x;
}
void BounceCalculator::setMaxN(int x){
	maxN = x;
}
void BounceCalculator::verboseOn(){
	verbose = true;
}
void BounceCalculator::verboseOff(){
	verbose = false;
}

}
