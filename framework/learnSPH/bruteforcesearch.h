
#pragma once
#include <vector>
#include <cmath>
#include <cstddef>
#include <Eigen/Dense>

namespace learnSPH
{
	namespace bfsearch
	{
		inline double distance(const Eigen::Vector3d &xi, const Eigen::Vector3d &xj)
		{
			const Eigen::Vector3d r = xi - xj;
			return r.norm();
		}
		inline void find_neighbors_short(const std::vector<Eigen::Vector3d> &particles,
								   double radius,
								   std::vector<std::vector<unsigned int> > &neighbors)
		{
			for(size_t i = 0; i<particles.size(); i++)
			{
				for(size_t j = 0; j<particles.size(); j++)
				{
					if (i == j) continue;
					//neighbors[i][j]= 
					if (distance(particles[i],particles[j]) <= radius) {
						neighbors[i].push_back(static_cast<unsigned int>(j));
					}
				}	
			}
		}

		inline void find_neighbors(const std::vector<Eigen::Vector3d> &particles,
								   double radius,
								   std::vector<std::vector<unsigned int> > &neighbors)
		{
			const size_t n = particles.size();
			neighbors.assign(n, std::vector<unsigned int>());
			if (n == 0) return;
			const double r2 = radius * radius;
			for (size_t i = 0; i < n; ++i) {
				neighbors[i].clear();
				const double xi = particles[i][0];
				const double yi = particles[i][1];
				const double zi = particles[i][2];
				for (size_t j = 0; j < n; ++j) {
					if (i == j) continue;
					const double dx = xi - particles[j][0];
					const double dy = yi - particles[j][1];
					const double dz = zi - particles[j][2];
					const double d2 = dx * dx + dy * dy + dz * dz;
					if (d2 <= r2) {
						neighbors[i].push_back(static_cast<unsigned int>(j));
					}
				}
			}
		}


	}
}
