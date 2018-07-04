#include "mfac.h"

#include <vector>
#include <cmath>
using std::vector;
//using std::cout; using std::endl;

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
	vector<double> y{ -1, 1, 0.5 };   //!< 实际输出，并初始化y(1), y(2), y(3)
	vector<double> u{ 0, 0 };         //!< 控制器输入
	vector<double> du{ 0, 0 };        //!< 控制器输入导数

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

	vector<double> phi{ 2, 2 };       //!< φ
	vector<double> a{ 0, 0 };
	vector<double> err{ 0, 0 };

	Mfac mfac(eta, mu, rho, lambda, 1e-5);
	for (int k = 2; k != 1000; ++k)
	{
		a.push_back(1 + round((k + 1) / 500.0));

		//使用控制器计算输入值u
		u.push_back(mfac.out(yd[k + 1], y[k]));

		// 根据模型计算输出值，实际中不需要
		if (k <= 499)
		{
			y.push_back(y[k] / (1 + pow(y[k], 2)) + pow(u[k], 3));
		}
		else
		{
			y.push_back((y[k] * y[k - 1] * y[k - 2] * u[k - 1] * (y[k - 2] - 1) + a[k] * u[k]) / (1 + pow(y[k - 1], 2) + pow(y[k - 2], 2)));
		}
		du.push_back(u[k] - u[k - 1]);

		err.push_back(yd[k] - y[k]);
	}

	return 0;
}