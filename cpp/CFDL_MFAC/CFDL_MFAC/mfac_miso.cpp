#include "mfac_miso.h"

/**
 * @brief Construct a new Mfac Miso:: Mfac Miso object
 * 
 * @note mfac的参数在这里设定，这不是一个有适用性的控制器，调试阶段为了好调试这么写。之后将参数设置为构造函数的变量，可以构造不同的控制器
 */
MfacMiso::MfacMiso() :
    eta_(1),
    mu_(1), 
    rho_(1),
    lambda_(3),
    epsilon_(1e-5),
    M_(50)
{
    y_.resize(2);
    u1_.resize(2);
    u2_.resize(2);
    phi1_.resize(2);
    phi2_.resize(2);

    y_[0] = 1;
    y_[1] = 0.5;
    
    u1_[0] = 1;
    u1_[1] = 1;

    u2_[0] = 1;
    u2_[1] = 1;

    phi1_init_value_ = phi1_[0] = 0.5;
    phi2_init_value_ = phi2_[0] = -0.2;
}

MfacMiso::~MfacMiso()
{
}

/**
 * @brief 控制器输出计算函数
 * 
 * @param yd 期望输出
 * @param y 实际输出
 * @return vector<double> 返回大小为2，第一个为u1，第二个为u2
 */
vector<double> MfacMiso::out(double yd, double y)
{
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