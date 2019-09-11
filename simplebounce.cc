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


////////////////////////////////////////////////////////////////////////////////
// class Scalarfield
/*
 * scalar field(s) with O(d) symmetry in d-dimensional Euclidean space
 */

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

// return the value of scalar field phi_iphi at r_i
double Scalarfield::phi(const int i, const int iphi) const {
	return phi_[i*nphi_+iphi];
}

// set the value of scalar field phi_iphi to phi__
void Scalarfield::setPhi(const int i, const int iphi, const double phi__) {
	phi_[i*nphi_+iphi] = phi__;
}

// add phi__ to the value of scalar field phi_iphi
void Scalarfield::addToPhi(const int i, const int iphi, const double phi__) {
	phi_[i*nphi_+iphi] += phi__;
}

// return the address of phi_0 at r_i. 
// phi_0, phi_1, phi_2, ... can be obtained by x[0], x[1], x[2], ... after x = phivec(i)
double* Scalarfield::phivec(const int i) const {
	return &phi_[i*nphi_+0];
}

// radius : r_i = i * dr
double Scalarfield::r(const int i) const {
	return dr_*i;
}

// Laplacian in radial coordinate : \nabla^2 \phi = d^2 phi / dr^2 + (d-1)/r * dphi/dr
// note that rinv_[i] = 1/r(i). See Eq. 9 in the manual
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
	if(rmax__ > 0.){
		rmax_ = rmax__;
	} else {
		std::cerr << "!!! rmax should be positive value !!!"<< std::endl;
		std::cerr << "!!! rmax is set to 1. !!!"<< std::endl;
		rmax_ = 1.;
	}
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

// set the number of scalar field(s)
void Scalarfield::setNphi(const int nphi__){
	nphi_ = nphi__;
	delete[] phi_;
	phi_ = new double[n_*nphi_];
}

// return the number of the grid
int Scalarfield::n() const{
	return n_;
}

// return the number of the scalar field(s)
int Scalarfield::nphi() const{
	return nphi_;
}

// return the dimension of space
int Scalarfield::dim() const{
	return dim_;
}

// return the radius at the boundary
double Scalarfield::rmax() const{
	return rmax_;
}

// return the lattice spacing
double Scalarfield::dr() const{
	return dr_;
}

// return pow(r(i), dim-1)
double Scalarfield::r_dminusoneth(const int i) const{
	return r_dminusoneth_[i];
}


////////////////////////////////////////////////////////////////////////////////
// class GenericModel
/*
 * This class should be overloaded by user. e.g.,
 *
	class MyModel : public GenericModel{
	  public:
		MyModel(){
			setNphi(1);
		}
		// potential for scalar field(s)
		double vpot (const double* phi) const{
			return phi[0]*phi[0]/2. - phi[0]*phi[0]*phi[0]/3.;
		}
		// first derivative(s) of potential
		void calcDvdphi(const double *phi, double *dvdphi) const{
			dvdphi[0] = phi[0] - phi[0]*phi[0];
		}
	};
 *
 */

// set the number of the scalar field(s)
void GenericModel::setNphi(const int nphi__){
	nphi_ = nphi__;
}

// return the number of scalar field(s)
int GenericModel::nphi() const{
	return nphi_;
}


////////////////////////////////////////////////////////////////////////////////
// class BounceCalculator
/*
 * calculate the bounce action by a gradient flow equation which is proposed by 1907.02417.
 * The gradient flow equation solves Coleman-Glaser-Martin's reduced problem. See, https://doi.org/10.1007/BF01609421
 * The bounce solution can be obtained from scale transformation. See rBounce() and action().
 */

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
	setNphi(model_->nphi());
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

// kinetic energy of the configuration
// \int_0^\infty dr r^{d-1} \sum_i (-1/2) \phi_i \nabla^2\phi_i
double BounceCalculator::t() const {
	double integrand[n()-1];
	for(int i=0; i<n()-1; i++){
		integrand[i] = 0.;
		for(int iphi=0; iphi<nphi(); iphi++){
			integrand[i] += r_dminusoneth(i) * -0.5 * phi(i,iphi) * lap(i,iphi);
		}
	}
	return integral(integrand, dr(), n()-1);
}

// potential energy of the configuration
// \int_0^\infty dr r^{d-1} V(\phi)
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

	// \nabla^2 \phi_iphi at r = r_i
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
		double dvdphi[nphi()];
		model->calcDvdphi(phivec(i), dvdphi);
		for(int iphi=0; iphi<nphi(); iphi++){
			integrand1[i] += r_dminusoneth(i) * dvdphi[iphi] * laplacian[i][iphi];
			integrand2[i] += r_dminusoneth(i) * dvdphi[iphi] * dvdphi[iphi];
		}
	}
	integrand1[n()-1] = 0.;
	integrand2[n()-1] = 0.;

	// Eq. 9 of 1907.02417
	lambda = integral(integrand1,dr(),n()) / integral(integrand2,dr(),n());

	// RHS of Eq. 8 of 1907.02417
	// phi at boundary is fixed to phiFV and will not be updated.
	double RHS[n()][nphi()];
	for(int i=0; i<n()-1; i++){
		double dvdphi[nphi()];
		model->calcDvdphi(phivec(i), dvdphi);
		for(int iphi=0; iphi<nphi(); iphi++){
			RHS[i][iphi] = laplacian[i][iphi] - lambda*dvdphi[iphi];
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
	double dvdphi[nphi()];
	model->calcDvdphi(phivec(i), dvdphi);
	return lap(i,iphi) - lambda*dvdphi[iphi];
}

// RHS of EOM for the bounce solution at r = sqrt(lambda) * r_i
// The bounce configuration can be obtained from Eq. 15 of 1907.02417
double BounceCalculator::residualBounce(const int i, const int iphi) const{
	double dvdphi[nphi()];
	model->calcDvdphi(phivec(i), dvdphi);
	return lap(i,iphi)/lambda - dvdphi[iphi];
}

// Kinetic energy of the bounce
// The bounce configuration can be obtained from Eq. 15 of 1907.02417
double BounceCalculator::tBounce() const {
	double area = dim() * pow(M_PI,dim()/2.) / tgamma(dim()/2.+1.);
	return area * pow(lambda, dim()/2.-1.)*t();
}

// Potential energy of the bounce
// The bounce configuration can be obtained from Eq. 15 of 1907.02417
double BounceCalculator::vBounce() const {
	double area = dim() * pow(M_PI,dim()/2.) / tgamma(dim()/2.+1.);
	return area * pow(lambda, dim()/2.)*v();
}

// Euclidean action in d-dimensional space 
// The bounce configuration can be obtained from Eq. 15 of 1907.02417
double BounceCalculator::action() const{
	return tBounce() + vBounce();
}

// this value should be one for the bounce solution
double BounceCalculator::oneIfBounce() const{
	return (2.-dim())/dim() * tBounce() / vBounce();
}

// boucne solution from scale transformation
double BounceCalculator::rBounce(const int i) const{
	return sqrt(lambda)*dr()*i;
}


// set the posiiton of a point which gives V < V_FV and false vacua
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
// See Eq. 11 in the manual.
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

// evolve the configuration from tau = 0 to tau = tauend
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
// See Fig. 1 of the manual
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

	bool finished = false;
	while(!finished){
		for(int i=0; i<tend1/tend0; i++){
			finished = true;
			double deriv = evolveUntil(tend0*rmax()*rmax());
			if(verbose){
				std::cerr << "\t" << "deriv :\t" << deriv << std::endl;
				std::cerr << "\t" << "field excursion :\t" << fieldExcursion() << std::endl;
				std::cerr << "\t" << "derivative at boundary:\t" << derivativeAtBoundary() << std::endl;
			}
			// if dphi/dr at the boundary is NOT small enough
			if( deriv > derivMax/rmax() ) {
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
				if(n()>maxN){
					std::cerr << "!!! n became too large !!!" << std::endl;
					return -1;
				}

				finished = false;
				break;
			}
		}
	}

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

// print the result
int BounceCalculator::printBounceDetails() const{
	if(!setModelDone){
		std::cerr << "!!! model has not been set yet !!!"<< std::endl;
		return -1;
	}
	if(!setVacuumDone){
		std::cerr << "!!! correct vacua have not been set yet !!!"<< std::endl;
		return -1;
	}

	std::cout << "Bounce: " << std::endl;
	std::cout << "\tS =\t" << action() << std::endl;
	std::cout << "\tT =\t" << tBounce() << std::endl;
	std::cout << "\tV =\t" << vBounce() << std::endl;

	std::cout << "\tr = 0 : ";
	std::cout << "(";
	for(int iphi=0; iphi<nphi()-1; iphi++){
		std::cout << phi(0,iphi) << ", ";
	}
	std::cout << phi(0,nphi()-1)<< ")\t";
	std::cout << "V = " << model->vpot(phivec(0)) << std::endl;

	std::cout << "\tr = " << rBounce(n()-1) << " : ";
	std::cout << "(";
	for(int iphi=0; iphi<nphi()-1; iphi++){
		std::cerr << phi(n()-1,iphi) << ", ";
	}
	std::cout << phi(n()-1,nphi()-1) << ")\t";
	std::cout << "V = " << model->vpot(phivec(n()-1)) << std::endl;

	std::cout << std::endl;

	std::cout << "Before rescaling: " << std::endl;
	std::cout << "\tT =\t" << t() << std::endl;
	std::cout << "\tV =\t" << v() << std::endl;
	std::cout << "\tlambda =\t" << lambda << std::endl;
	std::cout << std::endl;


	std::cout << "Goodness of the solution: " << std::endl;
	std::cout << "\tderiv\t" << derivativeAtBoundary()/fieldExcursion()*rmax() << std::endl;
	std::cout << "\t(2-d)/d*T/V =\t" << oneIfBounce() << std::endl;


	return 0;
}

void BounceCalculator::setSafetyfactor(double x){
	safetyfactor = x;
}

// set the parameter to determine the maximal variation of the field in evolve()
void BounceCalculator::setMaximumvariation(double x){
	maximumvariation = x;
}

// set the initial value of the parameter to determine the initail configuration
void BounceCalculator::setXTV0(double x){
	xTV0 = x;
}

// set the initial value of the parameter to determine the initail configuration
void BounceCalculator::setWidth0(double x){
	width0 = x;
}

// set the maximal value of the derivative of field at the boundary
void BounceCalculator::setDerivMax(double x){
	derivMax = x;
}

// set tau0
void BounceCalculator::setTend0(double x){
	tend0 = x;
}

// set tau1
void BounceCalculator::setTend1(double x){
	tend1 = x;
}

// set the maximal value of the grid
void BounceCalculator::setMaxN(int x){
	maxN = x;
}

// turn on verbose mode
void BounceCalculator::verboseOn(){
	verbose = true;
}

// turn off verbose mode
void BounceCalculator::verboseOff(){
	verbose = false;
}

}
