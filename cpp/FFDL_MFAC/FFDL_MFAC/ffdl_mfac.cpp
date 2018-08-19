#include "ffdl_mfac.h"
#include <numeric>
using std::accumulate;

int sign(double x)
{
	if (x> 0) return 1;
	if (x <0) return -1;
	return 0;
}

//double accumulate

const int Ly = 1;
const int Lu = 2;
double y_init[Ly + 1] = { 0, 0 };
double u_init[Lu + 1] = { 0, 0, 0 };
double phi_init[Ly + Lu] = { 0, 0, 0 };

/**
 * @brief Construct a new Mfac Miso:: Mfac Miso object
 * 
 * @note mfac的参数在这里设定，这不是一个有适用性的控制器，调试阶段为了好调试这么写。之后将参数设置为构造函数的变量，可以构造不同的控制器
 *       这里不考虑Ly <= 0 的情况
 */
FfdlMfac::FfdlMfac() :
    eta_(1),
    mu_(1), 
    lambda_(3),
    epsilon_(1e-5),
	Ly_(Ly),
	Lu_(Lu),
	phi_number_(Ly_ + Lu_)
{
	// 根据Ly和Lu定义需要记录的y和u的个数
    y_.resize(Ly_ + 1);
	dy_.resize(Ly_);
    u_.resize(Lu_ + 1);
	du_.resize(Lu_);

	// 初始化y、dy、u、du
	for (auto i = 0; i != y_.size(); ++i)
	{
		y_[i] = y_init[i];
	}
	for (auto i = 0; i != dy_.size(); ++i)
	{
		dy_[i] = y_[i + 1] - y_[i];
	}
	
	for (auto i = 0; i != u_.size(); ++i)
	{
		u_[i] = u_init[i];
	}
	for (auto i = 0; i != du_.size(); ++i)
	{
		du_[i] = u_[i + 1] - u_[i];
	}

	// 初始化phi
	phi_.resize(phi_number_);
	phi_init_value_.resize(phi_number_);
	for (auto i = 0; i != phi_.size(); ++i)
	{
		phi_[i] = phi_init[i];
		phi_init_value_[i] = phi_init[i];
	}
}

FfdlMfac::~FfdlMfac()
{
}

/**
 * @brief 控制器输出计算函数
 * 
 * @param yd 期望输出
 * @param y 实际输出
 * @return double 返回控制器输出
 */
double FfdlMfac::out(double yd, double y)
{
	vector<double> last_phi = phi_;		// 记录上一次的phi

	vector<double> h;					// h矩阵
	double h_2norm_2 = 0;				// h矩阵2范数的平方
	// 构造h矩阵
	for (auto yy : y_)
	{
		h.push_back(yy);
	}

	for (auto u : u_)
	{
		h.push_back(u);
	}

	for (auto hh : h)
	{
		h_2norm_2 += hh * hh;
	}

	// 计算新的phi
	double phi_h = 0;					// phi乘h矩阵的值
	for (auto i = 0; i != phi_number_; ++i)
	{
		phi_h += phi_[i] * h[i];
	}

	for (auto i = 0; i != phi_number_; ++i)
	{
		phi_[i] = last_phi[i] + eta_ * (y_[Ly_] - y_[Ly_ - 1] - phi_h) * h[i] / (mu_ + h_2norm_2);
	}

	// 判断是否重置phi
	double phi_2norm;					// phi的2范数
	for (auto p : phi_)
	{
		phi_2norm += p * p;
	}
	phi_2norm = sqrt(phi_2norm);

	if (phi_2norm < epsilon_ || sqrt(h_2norm_2) < epsilon_ || sign(phi_[Ly_ + 1]) != sign(phi_init_value_[Ly_ + 1]))
	{
		phi_ = phi_init_value_;
	}

	// 计算u
	u_[Lu_] = u_[Lu_ - 1] + (rho_[Ly_ + 1] * phi_[Ly_ + 1] * (yd - y_[Ly_]) - phi_[Ly_ + 1] * accumulate(phi_[0], phi_[Ly_ - 1], 0)

    double du1 = u1_[1] - u1_[0];
    double du2 = u2_[1] - u2_[0];

    y_[0] = y_[1];
    u1_[0] = u1_[1];
    u2_[0] = u2_[1];
    phi1_[0] = phi1_[1];
    phi2_[0] = phi2_[1];

    y_[1] = y;

    phi1_[1] = phi1_[0] + eta_ * (y_[1] - y_[0] - phi1_[0] * du1 - phi2_[0] * du2) * du1 / (mu_ + du1 * du1 + du2 * du2);
    phi2_[1] = phi2_[0] + eta_ * (y_[1] - y_[0] - phi1_[0] * du1 - phi2_[0] * du2) * du2 / (mu_ + du1 * du1 + du2 * du2);

    if (abs(phi1_[1]) < epsilon_ || abs(phi1_[1]) > M_)
    {
        phi1_[1] = phi1_init_value_;
    }
    if (abs(phi2_[1]) < epsilon_ || abs(phi2_[1]) > M_)
    {
        phi2_[1] = phi2_init_value_;
    }

    u1_[1] = u1_[0] + rho_ * phi1_[1] * (yd - y_[1]) / (lambda_ + phi1_[1] * phi1_[1] + phi2_[1] * phi2_[1]);
    u2_[1] = u2_[0] + rho_ * phi2_[1] * (yd - y_[1]) / (lambda_ + phi1_[1] * phi1_[1] + phi2_[1] * phi2_[1]);

    vector<double> ret_u;
    ret_u.push_back(u1_[1]); ret_u.push_back(u2_[1]);
    return ret_u;
}