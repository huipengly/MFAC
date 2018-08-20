/**
 * @brief 测试书上例题4.7
 * 
 * @file main.cpp
 * @author huipengly
 * @date 2018-08-19
 */
#include "ffdl_mfac.h"

#include <vector>
#include <cmath>
using std::vector;
//using std::cout; using std::endl;

int main()
{
    const int N = 701;
	double y[N]{ 0, 0, 0, 0, 1, 0.2, 0 }, yd[N]{ 0 };//, u[N]{0};
	double u[N] = { 0 };

	// 生成期望y
    for (int i = 1; i != N; ++i)
    {
		if (i <= 490)
		{
			yd[i] = 0.4 * pow(-1, round(i / 50.0));
		}
		else
		{
			yd[i] = 0.1 + 0.1 * pow(-1, round(i / 50.0));
		}
    }

	FfdlMfac mfac;

	// 仿真过程
	for (int i = 6; i != N - 1; ++i)
	{
		if (i == 8)
		{
			double a = 0;
		}
		u[i] = mfac.out(yd[i], y[i]);

		double a = 1, b = 0;
		b = 4 * round(i / 100.0) + sin(i / 100.0); 
		y[i + 1] = (-0.9*a*y[i] + (b + 1)*u[i]) / (1 + y[i] * y[i]);
	}

    return 0;
}