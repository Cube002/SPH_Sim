#include <Eigen/Dense>

#include "assignment5_scenes.h"
#include "simulation_runner.h"

int main(int argc, char **argv)
{
	using namespace learnsph_app;

	// Assignment 5 b: dam break in a rectangular tank.
	const SimulationSizes sizes = SimulationSizes::from_radius(0.005);
	const double enlargement = 3.5 * sizes.particle_radius;

	SimulationConfig config;
	config.name = "Assignment 5 b: dam break_LargeBoundary";
	config.output_prefix = "assignment5_b_dam_break_LargeBoundary";
	config.stiffness = 1000.0;
	config.fluid_viscosity = 0.0025;
	config.boundary_viscosity = 0.0;
	config.gravity = Eigen::Vector3d(0.0, 0.0, -9.81);
	config.max_time_step = 0.00025;
	config.simulation_time = 2.0;
	config.frame_interval = 0.01;
	config.fluid_mass_scale = 1.0;
	config.cull_outside_bounds = true;
	config.export_boundary_particles = false; //exportieren der boundarys zum debuggen.
	config.bounds_min = Eigen::Vector3d(
		-enlargement, -enlargement, -enlargement);
	config.bounds_max = Eigen::Vector3d(
		0.15 + enlargement, 0.8 + enlargement, 1.0 + enlargement);
	apply_command_line_overrides(config, argc, argv);

	run_wcsph_simulation(sizes, dam_break_scene(sizes), config);
	return 0;
}
