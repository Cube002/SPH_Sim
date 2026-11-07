#pragma once
#include <vector>

#include <Eigen/Dense>
#include <CompactNSearch/CompactNSearch>

namespace learnSPH
{
	namespace sph_math
	{
		void compute_density(
			const std::vector<Eigen::Vector3d> &particles,
			std::vector<double> &densities,
			const CompactNSearch::NeighborhoodSearch &nsearch,
			unsigned int point_set_id,
			double particle_mass,
			double h);
	};
};
