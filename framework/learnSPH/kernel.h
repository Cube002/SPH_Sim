#pragma once
#include <cassert>
#include <cmath>
#include <Eigen/Dense>

namespace learnSPH
{
	namespace kernel
	{
		constexpr double PI = 3.14159265358979323846;
		
		double cubic_spline(const double q);
		double cubic_grad_spline(const double q);
		
		// Kernel value W(xi - xj; h) = (1/h^3) * s(q), q = |xi-xj|/h
		inline double cubic_spline(const Eigen::Vector3d &xi, const Eigen::Vector3d &xj, const double h)
		{
			const Eigen::Vector3d r = xi - xj;
			const double dist = r.norm();
			if (h <= 0.0) return 0.0;
			const double q = dist / h;
			return cubic_spline(q) / (h * h * h);
		}
		
		// Gradient of kernel with respect to x_i: ∇_{xi} W(xi-xj; h)
		inline Eigen::Vector3d cubic_spline_grad(const Eigen::Vector3d &xi, const Eigen::Vector3d &xj, const double h)
		{
			const Eigen::Vector3d r = xi - xj;
			const double dist = r.norm();
			if (h <= 0.0) return Eigen::Vector3d::Zero();
			if (dist < 1e-18) return Eigen::Vector3d::Zero();
			const double q = dist / h;
			if (q >= 2.0) return Eigen::Vector3d::Zero();
			const double grad_q = cubic_grad_spline(q);
			// ∇W = (1/h^3) * s'(q) * ∇q = (1/h^3) * s'(q) * (r/(h*dist)) = s'(q)/(h^4 * dist) * r
			return (grad_q / (h * h * h * h * dist)) * r;
		}
	};
};