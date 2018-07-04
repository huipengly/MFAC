#pragma once
#include <cmath>

class Mfac
{
public:
	Mfac(const double & eta, const double & mu, const double & rho, const double & lambda, const double & epsilon);
	~Mfac();
	double out(double yd, double y);
protected:
	const double eta;           //!< η
	const double mu;            //!< μ
	const double rho;			//!< ρ
	const double lambda;		//!< λ
	const double epsilon;		//!< ε

	double yd;					//!< 期望输出值，书里的yd(k + 1)
	double y;					//!< 实际输出值
	double y_pre;				//!< 上一次的实际输出值
	double u;					//!< 控制器输入
	double u_pre;				//!< 上一次控制器输入
	double du;					//!< u的微分
	double du_pre;				//!< 上一次的u的微分
	double phi;					//!< PPD φ的估算值
	double phi_pre;				//!< 上一次PPD φ的估算值
	const double phi_first;		//!< φ(1)的值
};

