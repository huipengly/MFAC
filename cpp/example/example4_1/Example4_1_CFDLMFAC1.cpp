/**
* Copyright (C), huipengly, All rights reserved.
* @file Example4_1_CFDLMFAC1.cpp
* @brief C++实现无模型自适应控制书例题4.1
* @author huipengly
* @email huipengly@gmail.com
* @version 1.0
* @date 2018-07-03
* @update 
* @note 
* @warning 
* @todo 
* @history: 
*/

#include <iostream>
#include <vector>
#include <cmath>
using std::vector;
using std::cout; using std::endl;

int main()
{
    const double nu = 1;            //!< ν
    const double eta = 1;           //!< η
    const double mu = 1;            //!< μ
    const double rho = 0.6;         //!< ρ
    const double lambda = 0.1;      //!< λ
    const double pi = 3.14159;      //!< Π

    const int N = 1000;

    vector<double> yd;              //!< 期望输出
    vector<double> y{-1, 1, 0.5};   //!< 实际输出，并初始化y(1), y(2), y(3)
    vector<double> u{0, 0};         //!< 控制器输入
    vector<double> du{0, 0};        //!< 控制器输入导数

    //!计算期望值，实际中不需要
    for (int i = 0; i != 1001; ++i)
    {
        if (i <= 299)
        {
            yd.push_back(0.5 * pow((-1), round((i + 1) / 100.0)));
        }
        else if (i > 299 && i <= 699)
        {
            yd.push_back(0.5 * sin((i + 1) * pi / 100.0) + 0.3 * cos((i + 1) * pi / 50.0));
        }
        else if (i > 699)
        {
            yd.push_back(0.5 * pow((-1), round((i + 1) / 100.0)));
        }
    }

    vector<double> phi{2, 2};       //!< φ
    vector<double> a{0, 0};
    vector<double> err{0, 0};
    for (int i = 2; i != 1000; ++i)
    {
		a.push_back(1 + round((i + 1) / 500.0));
        phi.push_back(phi[i - 1] + eta * (y[i] - y[i - 1] - phi[i - 1] * du[i - 1]) * du[i - 1] / (mu + du[i - 1] * du[i - 1]));
        
        if (phi[i] < 1e-5 || sqrt(du[i - 1] * du[i - 1]) < 1e-5)
        {
            phi[i] = 0.5;
        }

        if (nu == 1)
        {
            u.push_back(u[i - 1] + rho * phi[i] * (yd[i + 1] - y[i]) / (lambda + phi[i] * phi[i]));
        }
        else
        {
			// 例4.1不需要，还没有实现
            // @todo 矩阵运算
            // u(k) = u(k-1)+rou*fai(k,1)*(yd(k+1)-y(k)-fai(k,2:nu)*du(k-1,1:nu-1)')/(lamda+fai(k,1).^2); 
            // u.push_back(u[i - 1] + rho * phi[i] * (yd[i + 1] - y(k) - phi[i] * du[i]))
        }

        // 根据模型计算输出值，实际中不需要
        if (i <= 499)
        {
            y.push_back(y[i] / (1 + pow(y[i], 2)) + pow(u[i], 3));
        }
        else
        {
            y.push_back((y[i] * y[i - 1] * y[i - 2] * u[i - 1] * (y[i - 2] - 1) + a[i] * u[i]) / (1 + pow(y[i - 1], 2) + pow(y[i - 2], 2)));
        }
        du.push_back(u[i] - u[i - 1]);

        err.push_back(yd[i] - y[i]);
    }

    return 0;
}