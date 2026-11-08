#pragma once
#include <vector>

#include <Eigen/Dense>
#include <CompactNSearch/CompactNSearch>

namespace learnSPH
{
	namespace sph_math
	{
		// Assignment 3: fluid density without boundary particles.
		void compute_density(
			const std::vector<Eigen::Vector3d> &fluid_particles,
			std::vector<double> &densities,
			const CompactNSearch::NeighborhoodSearch &nsearch,
			unsigned int fluid_point_set_id,
			double fluid_particle_mass,
			double h);

		// Assignment 3, Eq. (8): fluid density including boundary particles.
		void compute_density(
			const std::vector<Eigen::Vector3d> &fluid_particles,
			const std::vector<Eigen::Vector3d> &boundary_particles,
			std::vector<double> &densities,
			const CompactNSearch::NeighborhoodSearch &nsearch,
			unsigned int fluid_point_set_id,
			unsigned int boundary_point_set_id,
			double fluid_particle_mass,
			const std::vector<double> &boundary_volumes,
			double fluid_density,
			double h);

		// Assignment 2, Eq. (6): corrected volume of each boundary particle.
		void compute_boundary_volumes(
			const std::vector<Eigen::Vector3d> &boundary_particles,
			std::vector<double> &volumes,
			const CompactNSearch::NeighborhoodSearch &nsearch,
			unsigned int boundary_point_set_id,
			double h);

		// Assignment 2, Eq. (7): convert corrected volumes to boundary masses.
		void compute_boundary_masses(
			const std::vector<double> &volumes,
			std::vector<double> &masses,
			double boundary_density);

		// Assignment 2: compensated mass from a previously computed local density.
		void compute_compensated_mass(
			const std::vector<double> &densities,
			std::vector<double> &masses,
			double base_particle_mass,
			double boundary_density);

		// Assignment 4, Eq. (11): equation-of-state pressure, clamped to >= 0.
		void compute_pressure(
			const std::vector<double> &densities,
			std::vector<double> &pressures,
			double rest_density,
			double stiffness);

		// Assignment 4, Eq. (12): pressure acceleration without boundaries.
		// `accelerations` must be zeroed by the caller; this function adds its
		// contribution.
		void compute_pressure_acceleration(
			const std::vector<Eigen::Vector3d> &fluid_particles,
			const std::vector<double> &fluid_densities,
			const std::vector<double> &pressures,
			std::vector<Eigen::Vector3d> &accelerations,
			const CompactNSearch::NeighborhoodSearch &nsearch,
			unsigned int fluid_point_set_id,
			double fluid_particle_mass,
			double h);

		// Assignment 4, Eq. (12): pressure acceleration with fluid-boundary
		// interaction. `accelerations` must be zeroed by the caller.
		void compute_pressure_acceleration(
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
			double h);

		// Assignment 4, Eq. (13): friction acceleration without boundaries.
		// `accelerations` must be zeroed by the caller; this function adds its
		// contribution.
		void compute_friction_acceleration(
			const std::vector<Eigen::Vector3d> &fluid_particles,
			const std::vector<Eigen::Vector3d> &fluid_velocities,
			const std::vector<double> &fluid_densities,
			std::vector<Eigen::Vector3d> &accelerations,
			const CompactNSearch::NeighborhoodSearch &nsearch,
			unsigned int fluid_point_set_id,
			double fluid_particle_mass,
			double fluid_viscosity,
			double h);

		// Assignment 4, Eq. (13): friction acceleration with fluid-boundary
		// interaction. `boundary_velocities` may be empty for a stationary
		// boundary; if non-empty it must match `boundary_particles`. This already
		// supports future moving-boundary scenes. `accelerations` must be zeroed
		// by the caller.
		void compute_friction_acceleration(
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
			double h);

		// Assignment 4, Eq. (14): external/body acceleration.
		// `accelerations` must be zeroed by the caller; this function adds its
		// contribution.
		void compute_body_acceleration(
			std::vector<Eigen::Vector3d> &accelerations,
			const Eigen::Vector3d &body_acceleration);
	};
};
