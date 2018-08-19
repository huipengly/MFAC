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
double y_init[Ly + 1] = { 1, 0.2 };
double u_init[Lu + 1] = { 0, 0, 0.5 };
double phi_init[Ly + Lu] = { -2, 0.5, 0.2 };

/**
 * @brief Construct a new Mfac Miso:: Mfac Miso object
 * 
 * @note mfac的参数在这里设定，这不是一个有适用性的控制器，调试阶段为了好调试这么写。之后将参数设置为构造函数的变量，可以构造不同的控制器
 *       这里不考虑Ly <= 0 的情况
 */
FfdlMfac::FfdlMfac() :
    eta_(0.2),
    mu_(1), 
    lambda_(7),
    epsilon_(1e-5),
	rho_(0.7),
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
	// 更新y
	for (auto i = 0; i != y_.size() - 1; ++i)
	{
		y_[i] = y_[i + 1];
	}
	y_.back() = y;
	// 更新u
	for (auto i = 0; i != u_.size() - 1; ++i)
	{
		u_[i] = u_[i + 1];
	}

	vector<double> last_phi = phi_;		// 记录上一次的phi

	vector<double> h;					// h矩阵
	double h_2norm_2 = 0;				// h矩阵2范数的平方
	// 构造h矩阵
	for (auto yy : dy_)
	{
		h.push_back(yy);
	}

	//for (auto u : du_)
	//{
	//	h.push_back(u);
	//}

	for (auto i = du_.size() - 1; i != -1; --i)
	{
		h.push_back(du_[i]);
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
	//double phi_2norm = 0;					// phi的2范数
	//for (auto p : phi_)
	//{
	//	phi_2norm += p * p;
	//}
	//phi_2norm = sqrt(phi_2norm);

	//if (phi_2norm < epsilon_ || sqrt(h_2norm_2) < epsilon_ || sign(phi_[Ly_ + 1]) != sign(phi_init_value_[Ly_ + 1]))
	//{
	//	phi_ = phi_init_value_;
	//}

	// 更新dy
	for (auto i = 0; i != dy_.size() - 1; ++i)
	{
		dy_[i] = dy_[i + 1];
	}
	dy_.back() = y_[Ly_] - y_[Ly_ - 1];

	// 计算u
	//u_[Lu_] = u_[Lu_ - 1] + (rho_[Ly_ + 1 - 1] * phi_[Ly_ + 1 - 1] * (yd - y_[Ly_]) - phi_[Ly_ + 1] * accumulate(phi_[0], phi_[Ly_ - 1], 0) * 
	u_[Lu_] = u_[Lu_ - 1] + rho_ * phi_[Ly_ + 1 - 1] * (yd - y_[Ly_] - phi_[1 - 1] * dy_[Ly_ - 1] - phi_[Ly_ + Lu_ - 1] * du_[Lu_ - 1]) / (lambda_ + phi_[Ly_ + 1 - 1] * phi_[Ly_ + 1 - 1]);

	// 更新du
	for (auto i = 0; i != du_.size() - 1; ++i)
	{
		du_[i] = du_[i + 1];
	}
	du_.back() = u_[Lu_] - u_[Lu_ - 1];

	return u_.back();

}