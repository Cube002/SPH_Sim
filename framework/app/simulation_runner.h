#pragma once

#include <algorithm>
#include <cassert>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

#include <Eigen/Dense>
#include <CompactNSearch/CompactNSearch>

#include "../learnSPH/io.h"
#include "../learnSPH/parallel.h"
#include "../learnSPH/sph_math.h"

namespace learnsph_app
{
	struct SimulationSizes
	{
		double particle_radius;
		double fluid_sampling_distance;
		double boundary_sampling_distance;
		double smoothing_length;
		double compact_support;

		static SimulationSizes from_radius(double particle_radius)
		{
			const double particle_diameter = 2.0 * particle_radius;
			const double fluid_sampling_distance = particle_diameter;
			const double smoothing_length = 1.2 * fluid_sampling_distance;

			return {
				particle_radius,
				fluid_sampling_distance,
				0.8 * fluid_sampling_distance,
				smoothing_length,
				2.0 * smoothing_length
			};
		}
	};

	struct SimulationConfig
	{
		std::string name = "WCSPH simulation";
		std::string output_prefix = "simulation";
		double rest_density = 1000.0;
		double stiffness = 0.0;
		double fluid_viscosity = 0.0;
		double boundary_viscosity = 0.0;
		Eigen::Vector3d gravity = Eigen::Vector3d::Zero();
		double max_time_step = 0.001;
		double simulation_time = 0.1;
		double frame_interval = 0.02;
		double fluid_mass_scale = 1.0;
		double cfl_factor = 0.4;
		bool use_cfl = true;
		bool use_openmp = true;
		int num_threads = 0; // 0 = auto, capped by parallel::default_thread_cap

		bool cull_outside_bounds = false;
		Eigen::Vector3d bounds_min = Eigen::Vector3d::Zero();
		Eigen::Vector3d bounds_max = Eigen::Vector3d::Zero();
		double max_wall_clock_seconds = 3600.0;
		bool export_boundary_particles = false;
	};

	inline void apply_command_line_overrides(
		SimulationConfig &config, int argc, char **argv)
	{
		std::vector<std::string> positional;
		for (int i = 1; i < argc; ++i)
		{
			const std::string argument = argv[i];
			if (argument == "--export-boundary")
			{
				config.export_boundary_particles = true;
			}
			else if (argument == "--no-export-boundary")
			{
				config.export_boundary_particles = false;
			}
			else
			{
				positional.push_back(argument);
			}
		}
		if (positional.size() > 0) config.simulation_time = std::atof(positional[0].c_str());
		if (positional.size() > 1) config.num_threads = std::atoi(positional[1].c_str());
		if (positional.size() > 2) config.frame_interval = std::atof(positional[2].c_str());
	}

	struct SimulationScene
	{
		std::vector<Eigen::Vector3d> fluid_particles;
		std::vector<Eigen::Vector3d> fluid_velocities;
		std::vector<Eigen::Vector3d> boundary_particles;
		std::vector<Eigen::Vector3d> boundary_velocities;
	};

	inline void run_wcsph_simulation(
		const SimulationSizes &sizes,
		SimulationScene scene,
		const SimulationConfig &config)
	{
		if (scene.fluid_particles.empty())
		{
			std::cerr << config.name << ": no fluid particles" << std::endl;
			return;
		}

		const auto total_start = std::chrono::steady_clock::now();
		learnSPH::parallel::configure(
			config.use_openmp, config.num_threads);

		if (scene.fluid_velocities.empty())
		{
			scene.fluid_velocities.assign(
				scene.fluid_particles.size(), Eigen::Vector3d::Zero());
		}
		assert(scene.fluid_particles.size() == scene.fluid_velocities.size());

		const bool has_boundary = !scene.boundary_particles.empty();
		if (has_boundary)
		{
			assert(scene.boundary_velocities.empty()
				|| scene.boundary_particles.size()
					== scene.boundary_velocities.size());
		}

		const double base_fluid_mass = config.rest_density
			* sizes.fluid_sampling_distance
			* sizes.fluid_sampling_distance
			* sizes.fluid_sampling_distance;
		const double fluid_particle_mass =
			base_fluid_mass * config.fluid_mass_scale;

		CompactNSearch::NeighborhoodSearch nsearch(sizes.compact_support);
		const unsigned int fluid_id = nsearch.add_point_set(
			scene.fluid_particles[0].data(),
			scene.fluid_particles.size(), true);

		unsigned int boundary_id = 0;
		if (has_boundary)
		{
			boundary_id = nsearch.add_point_set(
				scene.boundary_particles[0].data(),
				scene.boundary_particles.size(), true);
		}
		nsearch.find_neighbors(true);

		// Assignment 2/5: corrected boundary volume and mass are computed once.
		// Moving boundaries would require refreshing these arrays after a
		// re-sampling step.
		std::vector<double> boundary_volumes;
		std::vector<double> boundary_masses;
		if (has_boundary)
		{
			learnSPH::sph_math::compute_boundary_volumes(
				scene.boundary_particles, boundary_volumes, nsearch,
				boundary_id, sizes.smoothing_length);
			learnSPH::sph_math::compute_boundary_masses(
				boundary_volumes, boundary_masses, config.rest_density);
		}

		const auto cull_outside_particles = [&]() -> std::size_t
		{
			if (!config.cull_outside_bounds) return 0;

			const auto is_inside = [&](const Eigen::Vector3d &x)
			{
				return x.allFinite()
					&& (x.array() >= config.bounds_min.array()).all()
					&& (x.array() <= config.bounds_max.array()).all();
			};

			std::size_t kept = 0;
			for (std::size_t i = 0; i < scene.fluid_particles.size(); ++i)
			{
				if (!is_inside(scene.fluid_particles[i])) continue;

				if (kept != i)
				{
					scene.fluid_particles[kept] = scene.fluid_particles[i];
					scene.fluid_velocities[kept] = scene.fluid_velocities[i];
				}
				++kept;
			}

			const std::size_t removed = scene.fluid_particles.size() - kept;
			if (removed == 0) return 0;

			scene.fluid_particles.resize(kept);
			scene.fluid_velocities.resize(kept);
			if (kept > 0)
			{
				nsearch.resize_point_set(
					fluid_id, scene.fluid_particles[0].data(), kept);
			}
			return removed;
		};

		std::vector<double> densities;
		std::vector<double> pressures;
		std::vector<Eigen::Vector3d> accelerations(
			scene.fluid_particles.size(), Eigen::Vector3d::Zero());

		const auto compute_density = [&]()
		{
			if (has_boundary)
			{
				learnSPH::sph_math::compute_density(
					scene.fluid_particles, scene.boundary_particles,
					densities, nsearch, fluid_id, boundary_id,
					fluid_particle_mass, boundary_volumes,
					config.rest_density, sizes.smoothing_length);
			}
			else
			{
				learnSPH::sph_math::compute_density(
					scene.fluid_particles, densities, nsearch, fluid_id,
					fluid_particle_mass, sizes.smoothing_length);
			}
		};

		const auto compute_accelerations = [&]()
		{
			accelerations.assign(
				scene.fluid_particles.size(), Eigen::Vector3d::Zero());

			if (has_boundary)
			{
				learnSPH::sph_math::compute_pressure_acceleration(
					scene.fluid_particles, scene.boundary_particles,
					densities, pressures, boundary_volumes, accelerations,
					nsearch, fluid_id, boundary_id, fluid_particle_mass,
					config.rest_density, sizes.smoothing_length);

				learnSPH::sph_math::compute_friction_acceleration(
					scene.fluid_particles, scene.boundary_particles,
					scene.fluid_velocities, scene.boundary_velocities,
					densities, boundary_volumes, accelerations,
					nsearch, fluid_id, boundary_id, fluid_particle_mass,
					config.fluid_viscosity, config.boundary_viscosity,
					sizes.smoothing_length);
			}
			else
			{
				learnSPH::sph_math::compute_pressure_acceleration(
					scene.fluid_particles, densities, pressures,
					accelerations, nsearch, fluid_id,
					fluid_particle_mass, sizes.smoothing_length);

				learnSPH::sph_math::compute_friction_acceleration(
					scene.fluid_particles, scene.fluid_velocities,
					densities, accelerations, nsearch, fluid_id,
					fluid_particle_mass, config.fluid_viscosity,
					sizes.smoothing_length);
			}

			learnSPH::sph_math::compute_body_acceleration(
				accelerations, config.gravity);
		};

		int frame_index = 0;
		const auto write_frame = [&]()
		{
			const std::string filename = "../res/" + config.output_prefix
				+ "_" + std::to_string(frame_index) + ".vtk";
			learnSPH::write_particles_to_vtk(
				filename, scene.fluid_particles, densities,
				scene.fluid_velocities);

			if (config.export_boundary_particles && has_boundary)
			{
				const std::string boundary_filename = "../res/"
					+ config.output_prefix + "_boundary_"
					+ std::to_string(frame_index) + ".vtk";
				learnSPH::write_particles_to_vtk(
					boundary_filename, scene.boundary_particles,
					boundary_masses, "mass");
			}
			++frame_index;
		};

		const auto compute_time_step = [&]()
		{
			double dt = config.max_time_step;
			if (!config.use_cfl) return dt;

			double max_speed = 0.0;
#ifdef _OPENMP
#pragma omp parallel for reduction(max:max_speed) schedule(static) if(learnSPH::parallel::should_run(scene.fluid_velocities.size()))
#endif
			for (std::size_t i = 0; i < scene.fluid_velocities.size(); ++i)
			{
				max_speed = std::max(
					max_speed, scene.fluid_velocities[i].norm());
			}

			// Eq. (10) has B = dp/drho, so the numerical sound speed is
			// sqrt(B). This CFL step keeps information from crossing more
			// than a fraction of the smoothing length per step.
			const double sound_speed =
				std::sqrt(std::max(config.stiffness, 0.0));
			const double denominator = sound_speed + max_speed;
			if (denominator > 1e-12)
			{
				dt = std::min(dt, config.cfl_factor
					* sizes.smoothing_length / denominator);
			}
			return dt;
		};

		std::cout << config.name << std::endl;
		std::cout << "  fluid particles: " << scene.fluid_particles.size()
			<< std::endl;
		std::cout << "  boundary particles: "
			<< scene.boundary_particles.size() << std::endl;
		if (config.export_boundary_particles && has_boundary)
		{
			std::cout << "  boundary export: enabled" << std::endl;
		}
		std::cout << "  OpenMP: "
			<< (learnSPH::parallel::enabled() ? "enabled" : "disabled")
			<< ", threads = " << learnSPH::parallel::thread_count()
			<< std::endl;
		std::cout << "  B = " << config.stiffness
			<< ", nu_f = " << config.fluid_viscosity
			<< ", nu_b = " << config.boundary_viscosity << std::endl;
		std::cout << "  gravity = (" << config.gravity[0] << ", "
			<< config.gravity[1] << ", " << config.gravity[2]
			<< ") m/s^2" << std::endl;
		std::cout << "  max dt = " << config.max_time_step
			<< " s, t_end = " << config.simulation_time << " s" << std::endl;

		// Initial state. Density is needed for the exported scalar field.
		std::size_t culled_particles = cull_outside_particles();
		if (scene.fluid_particles.empty())
		{
			std::cerr << config.name
				<< ": all fluid particles are outside the configured bounds"
				<< std::endl;
			return;
		}
		if (culled_particles > 0) nsearch.find_neighbors(true);
		compute_density();
		write_frame();

		const auto loop_start = std::chrono::steady_clock::now();
		double time = 0.0;
		double next_frame_time = config.frame_interval;
		int time_step = 0;
		bool stopped_by_wall_clock = false;

		while (time < config.simulation_time - 1e-12)
		{
			if (time_step % 100 == 0)
			{
				const auto now = std::chrono::steady_clock::now();
				const double elapsed =
					std::chrono::duration<double>(now - loop_start).count();
				if (elapsed >= config.max_wall_clock_seconds)
				{
					stopped_by_wall_clock = true;
					break;
				}
			}

			culled_particles += cull_outside_particles();
			if (scene.fluid_particles.empty())
			{
				std::cerr << config.name
					<< ": all fluid particles have been culled" << std::endl;
				break;
			}

			// Assignment 5, Algorithm 1: neighbors first, then density and forces.
			nsearch.find_neighbors(true);
			compute_density();
			learnSPH::sph_math::compute_pressure(
				densities, pressures, config.rest_density,
				config.stiffness);
			compute_accelerations();

			double dt = compute_time_step();
			dt = std::min(dt, config.simulation_time - time);
			if (dt <= 0.0) break;

			// Semi-implicit Euler integration.
#ifdef _OPENMP
#pragma omp parallel for schedule(static) if(learnSPH::parallel::should_run(scene.fluid_particles.size()))
#endif
			for (unsigned int i = 0; i < scene.fluid_particles.size(); ++i)
			{
				scene.fluid_velocities[i] += dt * accelerations[i];
				scene.fluid_particles[i] += dt * scene.fluid_velocities[i];
			}

			time += dt;
			++time_step;

			// Assignment 5b: velocity is included in every exported VTK frame.
			if (config.frame_interval > 0.0
				&& time + 1e-12 >= next_frame_time)
			{
				write_frame();
				next_frame_time += config.frame_interval;
			}
		}

		if (!scene.fluid_particles.empty())
		{
			double max_speed = 0.0;
			Eigen::Vector3d position_min =
				Eigen::Vector3d::Constant(std::numeric_limits<double>::infinity());
			Eigen::Vector3d position_max =
				Eigen::Vector3d::Constant(-std::numeric_limits<double>::infinity());
			for (std::size_t i = 0; i < scene.fluid_particles.size(); ++i)
			{
				max_speed = std::max(max_speed, scene.fluid_velocities[i].norm());
				position_min = position_min.cwiseMin(scene.fluid_particles[i]);
				position_max = position_max.cwiseMax(scene.fluid_particles[i]);
			}
			std::cout << "  final fluid particles: "
				<< scene.fluid_particles.size() << std::endl;
			std::cout << "  final max speed: " << max_speed << " m/s"
				<< std::endl;
			std::cout << "  final position bounds: ["
				<< position_min.transpose() << "] to ["
				<< position_max.transpose() << "]" << std::endl;
		}
		std::cout << "  culled particles: " << culled_particles << std::endl;
		if (stopped_by_wall_clock)
		{
			std::cout << "  stopped by wall-clock limit ("
				<< config.max_wall_clock_seconds << " s)" << std::endl;
		}

		const auto total_end = std::chrono::steady_clock::now();
		const double setup_seconds =
			std::chrono::duration<double>(loop_start - total_start).count();
		const double loop_seconds =
			std::chrono::duration<double>(total_end - loop_start).count();
		const double total_seconds =
			std::chrono::duration<double>(total_end - total_start).count();

		std::cout << "  finished at t = " << time << " s after "
			<< time_step << " steps" << std::endl;
		std::cout << "  setup time: " << setup_seconds << " s" << std::endl;
		std::cout << "  simulation loop time: " << loop_seconds << " s"
			<< std::endl;
		std::cout << "  total wall-clock time: " << total_seconds << " s"
			<< std::endl;
		std::cout << "  wrote " << frame_index
			<< " VTK frames to `<build_folder>/res`." << std::endl;
	}
}
