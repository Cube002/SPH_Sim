#include "sph_math.h"
#include "kernel.h"

#include <cassert>

void learnSPH::sph_math::compute_density(
	const std::vector<Eigen::Vector3d> &fluid_particles,
	std::vector<double> &densities,
	const CompactNSearch::NeighborhoodSearch &nsearch,
	unsigned int fluid_point_set_id,
	double fluid_particle_mass,
	double h)
{
	const auto &ps = nsearch.point_set(fluid_point_set_id);
	const unsigned int n = static_cast<unsigned int>(fluid_particles.size());
	densities.assign(n, 0.0);

	for (unsigned int i = 0; i < n; ++i)
	{
		// CompactNSearch does not return particle i as its own neighbor.
		double density = fluid_particle_mass
			* kernel::cubic_spline(fluid_particles[i], fluid_particles[i], h);

		const unsigned int num_neighbors = ps.n_neighbors(fluid_point_set_id, i);
		for (unsigned int k = 0; k < num_neighbors; ++k)
		{
			const unsigned int j = ps.neighbor(fluid_point_set_id, i, k);
			if (j == i) continue;

			density += fluid_particle_mass
				* kernel::cubic_spline(fluid_particles[i], fluid_particles[j], h);
		}

		densities[i] = density;
	}
}

void learnSPH::sph_math::compute_density(
	const std::vector<Eigen::Vector3d> &fluid_particles,
	const std::vector<Eigen::Vector3d> &boundary_particles,
	std::vector<double> &densities,
	const CompactNSearch::NeighborhoodSearch &nsearch,
	unsigned int fluid_point_set_id,
	unsigned int boundary_point_set_id,
	double fluid_particle_mass,
	const std::vector<double> &boundary_volumes,
	double fluid_density,
	double h)
{
	assert(boundary_particles.size() == boundary_volumes.size());

	const auto &ps = nsearch.point_set(fluid_point_set_id);
	const unsigned int n = static_cast<unsigned int>(fluid_particles.size());
	densities.assign(n, 0.0);

	for (unsigned int i = 0; i < n; ++i)
	{
		// Fluid contribution: m_f * sum_j W_ij, including the self term.
		double density = fluid_particle_mass
			* kernel::cubic_spline(fluid_particles[i], fluid_particles[i], h);

		const unsigned int num_fluid_neighbors =
			ps.n_neighbors(fluid_point_set_id, i);
		for (unsigned int k = 0; k < num_fluid_neighbors; ++k)
		{
			const unsigned int j = ps.neighbor(fluid_point_set_id, i, k);
			if (j == i) continue;

			density += fluid_particle_mass
				* kernel::cubic_spline(fluid_particles[i], fluid_particles[j], h);
		}

		// Boundary contribution: rho_f * sum_k V_k * W_ik.
		const unsigned int num_boundary_neighbors =
			ps.n_neighbors(boundary_point_set_id, i);
		for (unsigned int k = 0; k < num_boundary_neighbors; ++k)
		{
			const unsigned int boundary_index =
				ps.neighbor(boundary_point_set_id, i, k);

			density += fluid_density * boundary_volumes[boundary_index]
				* kernel::cubic_spline(
					fluid_particles[i], boundary_particles[boundary_index], h);
		}

		densities[i] = density;
	}
}

void learnSPH::sph_math::compute_boundary_volumes(
	const std::vector<Eigen::Vector3d> &boundary_particles,
	std::vector<double> &volumes,
	const CompactNSearch::NeighborhoodSearch &nsearch,
	unsigned int boundary_point_set_id,
	double h)
{
	const auto &ps = nsearch.point_set(boundary_point_set_id);
	const unsigned int n = static_cast<unsigned int>(boundary_particles.size());
	volumes.assign(n, 0.0);

	for (unsigned int i = 0; i < n; ++i)
	{
		// Assignment 2, Eq. (6): V_k = 1 / sum_l W_kl.
		double sum_kernel = kernel::cubic_spline(
			boundary_particles[i], boundary_particles[i], h);

		const unsigned int num_neighbors =
			ps.n_neighbors(boundary_point_set_id, i);
		for (unsigned int k = 0; k < num_neighbors; ++k)
		{
			const unsigned int j = ps.neighbor(boundary_point_set_id, i, k);
			if (j == i) continue;

			sum_kernel += kernel::cubic_spline(
				boundary_particles[i], boundary_particles[j], h);
		}

		volumes[i] = (sum_kernel > 0.0) ? (1.0 / sum_kernel) : 0.0;
	}
}

void learnSPH::sph_math::compute_boundary_masses(
	const std::vector<double> &volumes,
	std::vector<double> &masses,
	double boundary_density)
{
	masses.resize(volumes.size());

	// Assignment 2, Eq. (7): m_k = V_k * rho_b.
	for (std::size_t i = 0; i < volumes.size(); ++i)
	{
		masses[i] = volumes[i] * boundary_density;
	}
}

void learnSPH::sph_math::compute_compensated_mass(
	const std::vector<double> &densities,
	std::vector<double> &masses,
	double base_particle_mass,
	double boundary_density)
{
	masses.resize(densities.size());

	for (std::size_t i = 0; i < densities.size(); ++i)
	{
		// rho_k = m_b * sum_l W_kl  =>  V_k = m_b / rho_k.
		const double volume = (densities[i] > 0.0)
			? (base_particle_mass / densities[i])
			: 0.0;

		masses[i] = volume * boundary_density;
	}
}
void learnSPH::sph_math::compute_pressure_acceleration(
	const std::vector<Eigen::Vector3d> &fluid_particles,
	const std::vector<Eigen::Vector3d> &boundary_particles,
	const std::vector<double> &fluid_densities,
	const std::vector<double> &boundary_volumes,
	std::vector<double> &pressures,
	std::vector<Eigen::Vector3d> &accelerations,
	const CompactNSearch::NeighborhoodSearch &nsearch,
	unsigned int fluid_point_set_id,
	unsigned int boundary_point_set_id,
	double fluid_particle_mass,
	double rest_density,
	double stiffness,
	double h)
{}
void learnSPH::sph_math::compute_pressure_acceleration(
	const std::vector<Eigen::Vector3d> &fluid_particles,
	const std::vector<double> &fluid_densities,
	std::vector<Eigen::Vector3d> &accelerations,
	const CompactNSearch::NeighborhoodSearch &nsearch,
	unsigned int fluid_point_set_id,
	double fluid_particle_mass,
	double rest_density,
	double stiffness,
	double h)
{
	const auto &ps = nsearch.point_set(fluid_point_set_id);

	// HINT 1 -- Equation (11): resize `pressures` to fluid_particles.size()
	// and compute for every fluid particle i:
	//     pressures[i] = max(0.0, stiffness * (fluid_densities[i] - rest_density));
	int n = fluid_particles.size();
	std::vector<double> pressures(n,0.0);
	for(size_t i = 0; i < n; i++)//One hole loop just for computing pressures...
	{
		pressures[i]= std::max(0.0,stiffness*(fluid_densities[i] - rest_density));
	}

	for(size_t i = 0; i < n; i++)//n² loops
	{
		const unsigned int num_neighbors = ps.n_neighbors(fluid_point_set_id, i);
		Eigen::Vector3d acceleration(Eigen::Vector3d::Zero());
		for (unsigned int k = 0; k < num_neighbors; ++k)
		{
			const unsigned int j = ps.neighbor(fluid_point_set_id, i, k);
			if (j == i) continue;
			Eigen::Vector3d gradW_ij = kernel::cubic_grad_spline(fluid_particles[i], fluid_particles[j], h);
			acceleration += - fluid_particle_mass * (
				pressures[i]/(fluid_densities[i]*fluid_densities[i]) +
				pressures[j]/(fluid_densities[j]*fluid_densities[j]))
				* gradW_ij;			
		}
		accelerations[i] = acceleration;
	}
	// HINT 2 -- Equation (12), fluid-fluid part:
	// for each fluid neighbor j:
	//     a_i -= m_j * (p_i / rho_i^2 + p_j / rho_j^2) * gradW_ij
	// where gradW_ij = kernel::cubic_grad_spline(x_i, x_j, h).
	//
	// HINT 3 -- Equation (12), fluid-solid part:
	// for each boundary neighbor k:
	//     a_i -= rest_density * V_k * (p_i / rho_i^2) * gradW_ik
	// where V_k = boundary_volumes[k]. The boundary has no pressure of its own.
	//
	// HINT 4 -- Neighbor access:
	//     const auto &ps = nsearch.point_set(fluid_point_set_id);
	//     ps.n_neighbors(fluid_point_set_id, i)
	//     ps.neighbor(fluid_point_set_id, i, k)
	// Does CompactNSearch include particle i itself? Check the density code above.
	//
	// HINT 5 -- Add into accelerations[i], do not overwrite it. The pressure
	// kernel gradient at x_i - x_i is zero, so no explicit self term is needed.
}

void learnSPH::sph_math::compute_friction_acceleration(
	const std::vector<Eigen::Vector3d> &fluid_particles,
	const std::vector<Eigen::Vector3d> &boundary_particles,
	const std::vector<Eigen::Vector3d> &fluid_velocities,
	const std::vector<double> &fluid_densities,
	const std::vector<double> &boundary_volumes,
	std::vector<Eigen::Vector3d> &accelerations,
	const CompactNSearch::NeighborhoodSearch &nsearch,
	unsigned int fluid_point_set_id,
	unsigned int boundary_point_set_id,
	double fluid_particle_mass,
	double fluid_viscosity,
	double boundary_viscosity,
	double h)
{
	// HINT 1 -- Equation (13), fluid-fluid part:
	// for each fluid neighbor j:
	//     a_i += 2 * nu_f * (m_j / rho_j) * (v_i - v_j)
	//            * ((x_i - x_j).dot(gradW_ij))
	//            / (|x_i - x_j|^2 + 0.01 * h^2)
	// The dot product is a scalar; multiply it with the vector (v_i - v_j).
	//
	// HINT 2 -- Equation (13), fluid-solid part:
	// the boundary is stationary, so v_k = 0 and (v_i - v_k) = v_i:
	//     a_i += 2 * nu_b * V_k * v_i
	//            * ((x_i - x_k).dot(gradW_ik))
	//            / (|x_i - x_k|^2 + 0.01 * h^2)
	//
	// HINT 3 -- gradW_ij = kernel::cubic_grad_spline(x_i, x_j, h).
	// Use .squaredNorm() for |x_i - x_j|^2 and .dot() for the scalar product.
	//
	// HINT 4 -- Add into accelerations[i], do not overwrite it. The velocity
	// difference (v_i - v_i) is zero for j == i, so skipping self is harmless.
}

void learnSPH::sph_math::compute_body_acceleration(
	std::vector<Eigen::Vector3d> &accelerations,
	const Eigen::Vector3d &gravity)
{
	// HINT -- Equation (14): there are no neighbors involved here.
	// For every particle i: accelerations[i] += gravity;
	//
	// Typical gravity vector: Eigen::Vector3d(0.0, 0.0, -9.81).
	// This function should stay O(N) and needs no neighborhood search.
}
