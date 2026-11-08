#pragma once

#include <algorithm>
#include <cstddef>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace learnSPH
{
	namespace parallel
	{
		// Runtime switch shared by all SPH kernels. It is configured once by the
		// simulation runner before the time loop starts.
		inline bool g_enabled = true;
		inline constexpr std::size_t min_parallel_size = 128;
		inline constexpr int default_thread_cap = 8;

		inline void set_enabled(bool enabled)
		{
			g_enabled = enabled;
		}

		inline bool enabled()
		{
			return g_enabled;
		}

		inline bool should_run(std::size_t particle_count)
		{
			return g_enabled && particle_count >= min_parallel_size;
		}

		inline void configure(bool use_openmp, int num_threads)
		{
			set_enabled(use_openmp);
#ifdef _OPENMP
			if (!use_openmp)
			{
				omp_set_num_threads(1);
			}
			else
			{
				const int available_threads = omp_get_max_threads();
				const int threads = (num_threads > 0)
					? num_threads
					: std::min(available_threads, default_thread_cap);
				omp_set_num_threads(threads);
			}
#else
			(void)num_threads;
#endif
		}

		inline int thread_count()
		{
#ifdef _OPENMP
			return omp_get_max_threads();
#else
			return 1;
#endif
		}
	}
}
