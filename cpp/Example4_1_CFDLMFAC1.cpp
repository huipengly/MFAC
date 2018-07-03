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

int main()
{
    const double nu = 1;            //!<ν
    const double eta = 1;           //!<η
    const double mu = 1;            //!<μ
    const double rho = 0.6;         //!<ρ
    const double lambda = 0.1;      //!<λ
    const double pi = 3.14159;      //!<Π

    const int N = 1000;

    vector<double> yd;              //!<期望输出
    vector<double> y{-1, 1, 0.5};   //!<实际输出，并初始化y(1), y(2), y(3)
    vector<double> u{0, 0};         //!<控制器输入
    vector<double> du{0, 0};        //!<控制器输入导数

    //!计算期望值
    for (int i = 0; i != 1000; ++i)
    {
        if (i <= 300)
        {
            yd.push_back(0.5 * pow((-1), round(i / 100)));
        }
        else if (i > 300 && i <= 700)
        {
            yd.push_back(0.5 * sin(i * pi / 100) + 0.3 * cos(i * pi / 50));
        }
        else if (i > 700)
        {
            yd.push_back(0.5 * pow((-1), round( i / 100)));
        }
    }

    vector<double> phi{2, 2};       //!<φ
    vector<double> a{0, 0};
    for (int i = 2; i != 1000; ++i)
    {
        a.push_back(1 + round(i / 500));
        phi.push_back()
    }

    return 0;
}
