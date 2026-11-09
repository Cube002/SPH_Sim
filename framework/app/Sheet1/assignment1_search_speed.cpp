#include <stdlib.h>
#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <limits>
#include <cstddef>

#include <Eigen/Dense>
#include <CompactNSearch/CompactNSearch>

#include "../../learnSPH/io.h"
#include "../../learnSPH/sampling.h"
#include "../../learnSPH/bruteforcesearch.h"

using namespace std::chrono;

int main()
{
	std::cout << "Neighborhood Search Speed Test" << std::endl;

	// As per assignment: sample grid with particle distance = 2*r_particle
	// Use particle radius = 3mm, sampling_distance = 2*r = 6mm = 0.006 (or use 0.05 scaled; keep consistent)
	// We'll use sampling_distance = 0.006 for clarity, or adjust to get ~50x50x50
	const double r_particle = 0.003; // 3mm
	const double sampling_distance = 2.0 * r_particle; // particle distance = 2*r
	const double search_radius = 4.8 * r_particle; // compact support radius

	std::cout << "Sampling particles..." << std::endl;
	// Brute force test with smaller grid (50x50x50 is too big for O(n^2))
	std::vector<Eigen::Vector3d> particles_bf;
	{
		const Eigen::Vector3d bottom(0.0, 0.0, 0.0);
		const Eigen::Vector3d top(sampling_distance * 20.0, sampling_distance * 20.0, sampling_distance * 20.0);
		learnSPH::sampling::fluid_box(particles_bf, bottom, top, sampling_distance);
	}
	std::cout << "Particles (BF): " << particles_bf.size() << std::endl;

	// (a) Brute force
	{
		std::cout << "\n=== Brute Force ===" << std::endl;
		std::vector<std::vector<unsigned int> > neighbors_bf;
		// discard first run
		learnSPH::bfsearch::find_neighbors(particles_bf, search_radius, neighbors_bf);
		// timed runs
		double total_ms = 0.0;
		const int runs = 3;
		for (int run = 0; run < runs; ++run) {
			auto t0 = high_resolution_clock::now();
			learnSPH::bfsearch::find_neighbors(particles_bf, search_radius, neighbors_bf);
			auto t1 = high_resolution_clock::now();
			total_ms += duration<double, std::milli>(t1 - t0).count();
		}
		std::cout << "Avg time over " << runs << " runs: " << (total_ms / runs) << " ms" << std::endl;
		// average neighbors
		double sum_n = 0.0;
		for (const auto &nlist : neighbors_bf) sum_n += nlist.size();
		if (!neighbors_bf.empty()) std::cout << "Avg neighbors: " << (sum_n / neighbors_bf.size()) << std::endl;
	}

	// (a2) Brute force short
	{
		std::cout << "\n=== Brute Force Short ===" << std::endl;
		std::vector<std::vector<unsigned int> > neighbors_bfs(particles_bf.size());
		// discard first run
		learnSPH::bfsearch::find_neighbors_short(particles_bf, search_radius, neighbors_bfs);
		double total_ms = 0.0;
		const int runs = 3;
		for (int run = 0; run < runs; ++run) {
			for (auto &v : neighbors_bfs) v.clear();
			auto t0 = high_resolution_clock::now();
			learnSPH::bfsearch::find_neighbors_short(particles_bf, search_radius, neighbors_bfs);
			auto t1 = high_resolution_clock::now();
			total_ms += duration<double, std::milli>(t1 - t0).count();
		}
		std::cout << "Avg time over " << runs << " runs: " << (total_ms / runs) << " ms" << std::endl;
		double sum_n = 0.0;
		for (const auto &nlist : neighbors_bfs) sum_n += nlist.size();
		if (!neighbors_bfs.empty()) std::cout << "Avg neighbors: " << (sum_n / neighbors_bfs.size()) << std::endl;
	}

	// (b) CompactNSearch with 50x50x50 as reference
	{
		std::cout << "\n=== CompactNSearch ===" << std::endl;
		std::vector<Eigen::Vector3d> particles_cns;
		const Eigen::Vector3d bottom(0.0, 0.0, 0.0);
		const Eigen::Vector3d top(sampling_distance * 50.0, sampling_distance * 50.0, sampling_distance * 50.0);
		learnSPH::sampling::fluid_box(particles_cns, bottom, top, sampling_distance);
		std::cout << "Particles (CNS): " << particles_cns.size() << std::endl;
		CompactNSearch::NeighborhoodSearch nsearch(search_radius);
		// add point set (pass data pointer)
		if (!particles_cns.empty()) {
			nsearch.add_point_set(particles_cns[0].data(), particles_cns.size(), true);
			// discard first run
			nsearch.find_neighbors();
			double total_ms = 0.0;
			const int runs = 5;
			for (int run = 0; run < runs; ++run) {
				auto t0 = high_resolution_clock::now();
				nsearch.find_neighbors(true);
				auto t1 = high_resolution_clock::now();
				total_ms += duration<double, std::milli>(t1 - t0).count();
			}
			std::cout << "Avg time over " << runs << " runs: " << (total_ms / runs) << " ms" << std::endl;
			// average neighbors
			const auto &ps = nsearch.point_set(0);
			double sum_n = 0.0;
			for (size_t i = 0; i < ps.n_points(); ++i) {
				sum_n += ps.n_neighbors(0, static_cast<unsigned int>(i));
			}
			if (ps.n_points() > 0) std::cout << "Avg neighbors: " << (sum_n / ps.n_points()) << std::endl;
		}
	}

	return 0;
}