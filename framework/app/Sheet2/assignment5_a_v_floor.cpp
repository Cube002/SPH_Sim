#include <Eigen/Dense>

#include "assignment5_scenes.h"
#include "simulation_runner.h"

int main()
{
	using namespace learnsph_app;

	// Assignment 5 a(v): add a floor sampled with boundary particles.
	const SimulationSizes sizes = SimulationSizes::from_radius(0.03);

	SimulationConfig config;
	config.name = "Assignment 5 a(v): boundary floor";
	config.output_prefix = "assignment5_a_v_floor";
	config.stiffness = 1000.0;
	config.fluid_viscosity = 0.1;
	config.boundary_viscosity = 0.0;
	config.gravity = Eigen::Vector3d(0.0, 0.0, -9.81);
	config.max_time_step = 0.001;
	config.simulation_time = 0.4;
	config.frame_interval = 0.02;
	config.fluid_mass_scale = 1.0;

	run_wcsph_simulation(sizes, floor_scene(sizes), config);
	return 0;
}
