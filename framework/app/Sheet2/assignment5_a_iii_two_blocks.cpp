#include <Eigen/Dense>

#include "assignment5_scenes.h"
#include "simulation_runner.h"

int main()
{
	using namespace learnsph_app;

	// Assignment 5 a(iii): two 1 m^3 blocks collide at 2 m/s relative speed.
	const SimulationSizes sizes = SimulationSizes::from_radius(0.03);

	SimulationConfig config;
	config.name = "Assignment 5 a(iii): two colliding fluid blocks";
	config.output_prefix = "assignment5_a_iii_two_blocks";
	config.stiffness = 1000.0;
	config.fluid_viscosity = 0.1;
	config.boundary_viscosity = 0.0;
	config.max_time_step = 0.001;
	config.simulation_time = 1.0;
	config.frame_interval = 0.02;
	config.fluid_mass_scale = 0.9;

	run_wcsph_simulation(
		sizes, two_fluid_blocks_scene(sizes, 2.0), config);
	return 0;
}
