#include "mfac_miso.h"

#include <vector>
#include <cmath>
using std::vector;
//using std::cout; using std::endl;

int main()
{
    const int N = 1000;
    double y[1001]{1, 0.5, 0}, yd[1000]{0}, u1[1000]{0}, u2[1000]{0};

    vector<double> u(2);

    for (int i = 0; i != 1000; ++i)
    {
		yd[i] = pow(-1, round((i + 1) / 100.0));
    }

    MfacMiso mfac;

    for (int i = 1; i < 1000; ++i)
    {
		if (i == 49)
		{
			double aa = 1;
		}
        u = mfac.out(yd[i + 1], y[i]);
		u1[i] = u[0];
		u2[i] = u[1];
        y[i + 1] = (5 * y[i] + 2 * u1[i] - 3 * u2[i] * u2[i] + 2 * u1[i] * u1[i]) / (5 + u1[i] + 5 * u2[i]);
    }

    return 0;
}