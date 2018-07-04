#include "mfac.h"



Mfac::Mfac(const double &eta, const double &mu, const double &rho, const double &lambda, const double &epsilon) :
	eta(eta),
	mu(mu), 
	rho(rho),
	lambda(lambda),
	epsilon(epsilon),
	y_pre(1),
	u_pre(0),
	u(0),
	du(0),
	du_pre(0),
	phi_first(2),
	phi_pre(2)		//phi_pre(phi_first)为什么错，C++ prime上看到过
{
}

Mfac::~Mfac()
{
}

double Mfac::out(double yd, double y)
{
	//更新期望值和实际输出值
	// @todo:研究函数变量和成员函数同名时怎么办
	this->yd = yd;
	y_pre = this->y;
	this->y = y;

	phi = phi_pre + eta * du_pre  * (y - y_pre - phi_pre * du_pre)/ (mu + pow(du_pre, 2));

	if (fabs(phi) < epsilon || fabs(du_pre) < epsilon)
	{
		phi = phi_first;
	}

	u = u_pre + rho * phi * (yd - y) / (lambda + pow(phi, 2));

	phi_pre = phi;
	du_pre = u - u_pre;
	u_pre = u;

	return u;
}