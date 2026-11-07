#include "sph_math.h"
#include "kernel.h"

void learnSPH::sph_math::compute_density(
	const std::vector<Eigen::Vector3d> &particles,
	std::vector<double> &densities,
	const CompactNSearch::NeighborhoodSearch &nsearch,
	unsigned int point_set_id,
	double particle_mass,
	double h)
{
	const auto &ps = nsearch.point_set(point_set_id);
	const unsigned int n = static_cast<unsigned int>(particles.size());
	densities.assign(n, 0.0);

	for (unsigned int i = 0; i < n; ++i)
	{
		// CompactNSearch does not return particle i as its own neighbor.
		double density = particle_mass
			* kernel::cubic_spline(particles[i], particles[i], h);

		const unsigned int num_neighbors = ps.n_neighbors(point_set_id, i);
		for (unsigned int k = 0; k < num_neighbors; ++k)
		{
			const unsigned int j = ps.neighbor(point_set_id, i, k);
			if (j == i) continue;

			density += particle_mass
				* kernel::cubic_spline(particles[i], particles[j], h);
		}

		densities[i] = density;
	}
}
