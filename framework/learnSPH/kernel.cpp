#include "kernel.h"

double learnSPH::kernel::cubic_spline(const double q)
{
	assert(q >= 0.0);
	constexpr double alpha = 3.0 / (2.0 * PI);

	if (q < 1.0) {
		return alpha * (2.0/3.0 - q*q + 0.5*q*q*q);
	}
	else if (q < 2.0) {		
		return alpha * 1.0/6.0 * std::pow(2-q,3);
	}
	else {
		return 0;
	}
}	

double learnSPH::kernel::cubic_grad_spline(const double q)
{
	assert(q >= 0.0);
	constexpr double alpha = 3.0 / (2.0 * PI);

	if (q < 1.0) {
		// s'(q) = alpha * (-2q + 1.5*q^2)
		return alpha * (-2.0 * q + 1.5 * q * q);
	}
	else if (q < 2.0) {
		// s'(q) = d/dq [alpha/6*(2-q)^3] = alpha/6*3*(2-q)^2*(-1) = -alpha/2*(2-q)^2
		const double dq = 2.0 - q;
		return -alpha * 0.5 * dq * dq;
	}
	else {
		return 0.0;
	}
}

