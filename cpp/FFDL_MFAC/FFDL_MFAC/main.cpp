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

	for (int i = 6; i != N - 1; ++i)
	{
		if (i == 475)
		{
			double a = 0;
		}
		u[i] = mfac.out(yd[i], y[i]);

		double a = 1, b = 0;
		b = 4 * round(i / 100.0) + sin(i / 100.0); 
		y[i + 1] = (-0.9*a*y[i] + (b + 1)*u[i]) / (1 + y[i] * y[i]);
	}

  //  for (int i = 1; i < 1000; ++i)
  //  {
		//if (i == 49)
		//{
		//	double aa = 1;
		//}
  //      u = mfac.out(yd[i + 1], y[i]);
		//u1[i] = u[0];
		//u2[i] = u[1];
  //      y[i + 1] = (5 * y[i] + 2 * u1[i] - 3 * u2[i] * u2[i] + 2 * u1[i] * u1[i]) / (5 + u1[i] + 5 * u2[i]);
  //  }

    return 0;
}