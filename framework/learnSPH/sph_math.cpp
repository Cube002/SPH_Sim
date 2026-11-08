#include "sph_math.h"
#include "kernel.h"
#include "parallel.h"

#include <algorithm>
#include <cassert>

namespace
{
	void add_fluid_pressure_acceleration(
		unsigned int i,
		const std::vector<Eigen::Vector3d> &fluid_particles,
		const std::vector<double> &fluid_densities,
		const std::vector<double> &pressures,
		const CompactNSearch::NeighborhoodSearch &nsearch,
		unsigned int fluid_point_set_id,
		double fluid_particle_mass,
		double h,
		Eigen::Vector3d &acceleration)
	{
		const double density_i = fluid_densities[i];
		if (density_i <= 0.0) return;

		const auto &ps = nsearch.point_set(fluid_point_set_id);
		const double pressure_term_i =
			pressures[i] / (density_i * density_i);

		const unsigned int num_neighbors =
			ps.n_neighbors(fluid_point_set_id, i);
		for (unsigned int k = 0; k < num_neighbors; ++k)
		{
			const unsigned int j = ps.neighbor(fluid_point_set_id, i, k);
			if (j == i) continue;

			const double density_j = fluid_densities[j];
			if (density_j <= 0.0) continue;

			const double pressure_term_j =
				pressures[j] / (density_j * density_j);
			const Eigen::Vector3d grad_w = learnSPH::kernel::cubic_grad_spline(
				fluid_particles[i], fluid_particles[j], h);

			acceleration -= fluid_particle_mass
				* (pressure_term_i + pressure_term_j) * grad_w;
		}
	}

	void add_fluid_friction_acceleration(
		unsigned int i,
		const std::vector<Eigen::Vector3d> &fluid_particles,
		const std::vector<Eigen::Vector3d> &fluid_velocities,
		const std::vector<double> &fluid_densities,
		const CompactNSearch::NeighborhoodSearch &nsearch,
		unsigned int fluid_point_set_id,
		double fluid_particle_mass,
		double fluid_viscosity,
		double h,
		Eigen::Vector3d &acceleration)
	{
		const auto &ps = nsearch.point_set(fluid_point_set_id);
		const double denominator_offset = 0.01 * h * h;

		const unsigned int num_neighbors =
			ps.n_neighbors(fluid_point_set_id, i);
		for (unsigned int k = 0; k < num_neighbors; ++k)
		{
			const unsigned int j = ps.neighbor(fluid_point_set_id, i, k);
			if (j == i) continue;

			const double density_j = fluid_densities[j];
			if (density_j <= 0.0) continue;

			const Eigen::Vector3d displacement =
				fluid_particles[i] - fluid_particles[j];
			const double denominator =
				displacement.squaredNorm() + denominator_offset;
			if (denominator <= 0.0) continue;

			const Eigen::Vector3d grad_w = learnSPH::kernel::cubic_grad_spline(
				fluid_particles[i], fluid_particles[j], h);
			const double directional_term = displacement.dot(grad_w)
				/ denominator;

			acceleration += 2.0 * fluid_viscosity
				* (fluid_particle_mass / density_j)
				* directional_term
				* (fluid_velocities[i] - fluid_velocities[j]);
		}
	}
}

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

#ifdef _OPENMP
#pragma omp parallel for schedule(static) if(learnSPH::parallel::should_run(n))
#endif
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

#ifdef _OPENMP
#pragma omp parallel for schedule(static) if(learnSPH::parallel::should_run(n))
#endif
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

#ifdef _OPENMP
#pragma omp parallel for schedule(static) if(learnSPH::parallel::should_run(n))
#endif
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
#ifdef _OPENMP
#pragma omp parallel for schedule(static) if(learnSPH::parallel::should_run(volumes.size()))
#endif
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

#ifdef _OPENMP
#pragma omp parallel for schedule(static) if(learnSPH::parallel::should_run(densities.size()))
#endif
	for (std::size_t i = 0; i < densities.size(); ++i)
	{
		// rho_k = m_b * sum_l W_kl  =>  V_k = m_b / rho_k.
		const double volume = (densities[i] > 0.0)
			? (base_particle_mass / densities[i])
			: 0.0;

		masses[i] = volume * boundary_density;
	}
}

void learnSPH::sph_math::compute_pressure(
	const std::vector<double> &densities,
	std::vector<double> &pressures,
	double rest_density,
	double stiffness)
{
	pressures.resize(densities.size());

	// Assignment 4, Eq. (11): clamp negative pressures to zero.
#ifdef _OPENMP
#pragma omp parallel for schedule(static) if(learnSPH::parallel::should_run(densities.size()))
#endif
	for (std::size_t i = 0; i < densities.size(); ++i)
	{
		pressures[i] = std::max(
			0.0, stiffness * (densities[i] - rest_density));
	}
}

void learnSPH::sph_math::compute_pressure_acceleration(
	const std::vector<Eigen::Vector3d> &fluid_particles,
	const std::vector<double> &fluid_densities,
	const std::vector<double> &pressures,
	std::vector<Eigen::Vector3d> &accelerations,
	const CompactNSearch::NeighborhoodSearch &nsearch,
	unsigned int fluid_point_set_id,
	double fluid_particle_mass,
	double h)
{
	assert(fluid_particles.size() == fluid_densities.size());
	assert(fluid_particles.size() == pressures.size());
	if (accelerations.size() != fluid_particles.size())
	{
		accelerations.assign(
			fluid_particles.size(), Eigen::Vector3d::Zero());
	}

#ifdef _OPENMP
#pragma omp parallel for schedule(static) if(learnSPH::parallel::should_run(fluid_particles.size()))
#endif
	for (unsigned int i = 0; i < fluid_particles.size(); ++i)
	{
		Eigen::Vector3d acceleration = Eigen::Vector3d::Zero();
		add_fluid_pressure_acceleration(
			i, fluid_particles, fluid_densities, pressures,
			nsearch, fluid_point_set_id, fluid_particle_mass, h,
			acceleration);
		accelerations[i] += acceleration;
	}
}

void learnSPH::sph_math::compute_pressure_acceleration(
	const std::vector<Eigen::Vector3d> &fluid_particles,
	const std::vector<Eigen::Vector3d> &boundary_particles,
	const std::vector<double> &fluid_densities,
	const std::vector<double> &pressures,
	const std::vector<double> &boundary_volumes,
	std::vector<Eigen::Vector3d> &accelerations,
	const CompactNSearch::NeighborhoodSearch &nsearch,
	unsigned int fluid_point_set_id,
	unsigned int boundary_point_set_id,
	double fluid_particle_mass,
	double rest_density,
	double h)
{
	assert(fluid_particles.size() == fluid_densities.size());
	assert(fluid_particles.size() == pressures.size());
	assert(boundary_particles.size() == boundary_volumes.size());
	if (accelerations.size() != fluid_particles.size())
	{
		accelerations.assign(
			fluid_particles.size(), Eigen::Vector3d::Zero());
	}

	const auto &ps = nsearch.point_set(fluid_point_set_id);
#ifdef _OPENMP
#pragma omp parallel for schedule(static) if(learnSPH::parallel::should_run(fluid_particles.size()))
#endif
	for (unsigned int i = 0; i < fluid_particles.size(); ++i)
	{
		const double density_i = fluid_densities[i];
		if (density_i <= 0.0) continue;

		Eigen::Vector3d acceleration = Eigen::Vector3d::Zero();
		add_fluid_pressure_acceleration(
			i, fluid_particles, fluid_densities, pressures,
			nsearch, fluid_point_set_id, fluid_particle_mass, h,
			acceleration);

		const double pressure_term_i =
			pressures[i] / (density_i * density_i);
		const unsigned int num_boundary_neighbors =
			ps.n_neighbors(boundary_point_set_id, i);
		for (unsigned int k = 0; k < num_boundary_neighbors; ++k)
		{
			const unsigned int boundary_index =
				ps.neighbor(boundary_point_set_id, i, k);
			const Eigen::Vector3d grad_w = kernel::cubic_grad_spline(
				fluid_particles[i], boundary_particles[boundary_index], h);

			acceleration -= rest_density
				* boundary_volumes[boundary_index]
				* pressure_term_i * grad_w;
		}

		accelerations[i] += acceleration;
	}
}

void learnSPH::sph_math::compute_friction_acceleration(
	const std::vector<Eigen::Vector3d> &fluid_particles,
	const std::vector<Eigen::Vector3d> &fluid_velocities,
	const std::vector<double> &fluid_densities,
	std::vector<Eigen::Vector3d> &accelerations,
	const CompactNSearch::NeighborhoodSearch &nsearch,
	unsigned int fluid_point_set_id,
	double fluid_particle_mass,
	double fluid_viscosity,
	double h)
{
	assert(fluid_particles.size() == fluid_velocities.size());
	assert(fluid_particles.size() == fluid_densities.size());
	if (accelerations.size() != fluid_particles.size())
	{
		accelerations.assign(
			fluid_particles.size(), Eigen::Vector3d::Zero());
	}

#ifdef _OPENMP
#pragma omp parallel for schedule(static) if(learnSPH::parallel::should_run(fluid_particles.size()))
#endif
	for (unsigned int i = 0; i < fluid_particles.size(); ++i)
	{
		Eigen::Vector3d acceleration = Eigen::Vector3d::Zero();
		add_fluid_friction_acceleration(
			i, fluid_particles, fluid_velocities, fluid_densities,
			nsearch, fluid_point_set_id, fluid_particle_mass,
			fluid_viscosity, h, acceleration);
		accelerations[i] += acceleration;
	}
}

void learnSPH::sph_math::compute_friction_acceleration(
	const std::vector<Eigen::Vector3d> &fluid_particles,
	const std::vector<Eigen::Vector3d> &boundary_particles,
	const std::vector<Eigen::Vector3d> &fluid_velocities,
	const std::vector<Eigen::Vector3d> &boundary_velocities,
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
	assert(fluid_particles.size() == fluid_velocities.size());
	assert(fluid_particles.size() == fluid_densities.size());
	assert(boundary_particles.size() == boundary_volumes.size());
	assert(boundary_velocities.empty()
		|| boundary_particles.size() == boundary_velocities.size());
	if (accelerations.size() != fluid_particles.size())
	{
		accelerations.assign(
			fluid_particles.size(), Eigen::Vector3d::Zero());
	}

	const auto &ps = nsearch.point_set(fluid_point_set_id);
	const double denominator_offset = 0.01 * h * h;

#ifdef _OPENMP
#pragma omp parallel for schedule(static) if(learnSPH::parallel::should_run(fluid_particles.size()))
#endif
	for (unsigned int i = 0; i < fluid_particles.size(); ++i)
	{
		Eigen::Vector3d acceleration = Eigen::Vector3d::Zero();
		add_fluid_friction_acceleration(
			i, fluid_particles, fluid_velocities, fluid_densities,
			nsearch, fluid_point_set_id, fluid_particle_mass,
			fluid_viscosity, h, acceleration);

		const unsigned int num_boundary_neighbors =
			ps.n_neighbors(boundary_point_set_id, i);
		for (unsigned int k = 0; k < num_boundary_neighbors; ++k)
		{
			const unsigned int boundary_index =
				ps.neighbor(boundary_point_set_id, i, k);

			const Eigen::Vector3d displacement =
				fluid_particles[i] - boundary_particles[boundary_index];
			const double denominator =
				displacement.squaredNorm() + denominator_offset;
			if (denominator <= 0.0) continue;

			const Eigen::Vector3d grad_w = kernel::cubic_grad_spline(
				fluid_particles[i], boundary_particles[boundary_index], h);
			const double directional_term = displacement.dot(grad_w)
				/ denominator;

			const Eigen::Vector3d boundary_velocity =
				boundary_velocities.empty()
					? Eigen::Vector3d::Zero()
					: boundary_velocities[boundary_index];

			acceleration += 2.0 * boundary_viscosity
				* boundary_volumes[boundary_index]
				* directional_term
				* (fluid_velocities[i] - boundary_velocity);
		}

		accelerations[i] += acceleration;
	}
}

void learnSPH::sph_math::compute_body_acceleration(
	std::vector<Eigen::Vector3d> &accelerations,
	const Eigen::Vector3d &body_acceleration)
{
	// Assignment 4, Eq. (14): add the external/body acceleration.
#ifdef _OPENMP
#pragma omp parallel for schedule(static) if(learnSPH::parallel::should_run(accelerations.size()))
#endif
	for (std::size_t i = 0; i < accelerations.size(); ++i)
	{
		accelerations[i] += body_acceleration;
	}
}
