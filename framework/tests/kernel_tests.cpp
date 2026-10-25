#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include "../learnSPH/kernel.h"
#include <cmath>
#include <random>

namespace k = learnSPH::kernel;
constexpr double PI = 3.14159265358979323846;
constexpr double alpha = 3.0 / (2.0 * PI);

TEST_CASE("cubic_spline_grad: finite differences", "[kernel][grad]")
{
    const double h = 0.5;
    Eigen::Vector3d xi(0.1, 0.2, 0.3); //two arbitrary points so that the diff lies in the first branch
    Eigen::Vector3d xj(0.05, 0.15, 0.25);
    
    Eigen::Vector3d grad_ana = k::cubic_spline_grad(xi, xj, h); //The function to test
    
    const double eps = 1e-6;
    Eigen::Vector3d grad_fd = Eigen::Vector3d::Zero();
    
    // x component
    {
        Eigen::Vector3d xi1 = xi + Eigen::Vector3d(eps, 0, 0);
        Eigen::Vector3d xi2 = xi - Eigen::Vector3d(eps, 0, 0);
        grad_fd[0] = (k::cubic_spline(xi1, xj, h) - k::cubic_spline(xi2, xj, h)) / (2.0 * eps);//needs 3 components because the grad of kernel is direction dependant.
    }
    {
        Eigen::Vector3d xi1 = xi + Eigen::Vector3d(0, eps, 0);
        Eigen::Vector3d xi2 = xi - Eigen::Vector3d(0, eps, 0);
        grad_fd[1] = (k::cubic_spline(xi1, xj, h) - k::cubic_spline(xi2, xj, h)) / (2.0 * eps);
    }
    {
        Eigen::Vector3d xi1 = xi + Eigen::Vector3d(0, 0, eps);
        Eigen::Vector3d xi2 = xi - Eigen::Vector3d(0, 0, eps);
        grad_fd[2] = (k::cubic_spline(xi1, xj, h) - k::cubic_spline(xi2, xj, h)) / (2.0 * eps);
    }
    
    CHECK(grad_ana[0] == Catch::Approx(grad_fd[0]).epsilon(1e-6));
    CHECK(grad_ana[1] == Catch::Approx(grad_fd[1]).epsilon(1e-6));
    CHECK(grad_ana[2] == Catch::Approx(grad_fd[2]).epsilon(1e-6));
}

TEST_CASE("cubic_spline_grad: branch 2 finite differences", "[kernel][grad]")
{
    const double h = 0.3;
    Eigen::Vector3d xi(0.0, 0.0, 0.0);
    Eigen::Vector3d xj(-0.4, 0.2, 0.1); // q will be in [1,2)
    
    Eigen::Vector3d grad_ana = k::cubic_spline_grad(xi, xj, h);
    
    const double eps = 1e-6;
    Eigen::Vector3d grad_fd = Eigen::Vector3d::Zero();
    
    {
        Eigen::Vector3d xi1 = xi + Eigen::Vector3d(eps, 0, 0);
        Eigen::Vector3d xi2 = xi - Eigen::Vector3d(eps, 0, 0);
        grad_fd[0] = (k::cubic_spline(xi1, xj, h) - k::cubic_spline(xi2, xj, h)) / (2.0 * eps);
    }
    {
        Eigen::Vector3d xi1 = xi + Eigen::Vector3d(0, eps, 0);
        Eigen::Vector3d xi2 = xi - Eigen::Vector3d(0, eps, 0);
        grad_fd[1] = (k::cubic_spline(xi1, xj, h) - k::cubic_spline(xi2, xj, h)) / (2.0 * eps);
    }
    {
        Eigen::Vector3d xi1 = xi + Eigen::Vector3d(0, 0, eps);
        Eigen::Vector3d xi2 = xi - Eigen::Vector3d(0, 0, eps);
        grad_fd[2] = (k::cubic_spline(xi1, xj, h) - k::cubic_spline(xi2, xj, h)) / (2.0 * eps);
    }
    
    CHECK(grad_ana[0] == Catch::Approx(grad_fd[0]).epsilon(1e-6));
    CHECK(grad_ana[1] == Catch::Approx(grad_fd[1]).epsilon(1e-6));
    CHECK(grad_ana[2] == Catch::Approx(grad_fd[2]).epsilon(1e-6));
}

TEST_CASE("cubic_spline_grad: zero distance", "[kernel][grad]")
{
    Eigen::Vector3d xi(1.0, 2.0, 3.0);
    Eigen::Vector3d xj(1.0, 2.0, 3.0);
    const double h = 0.5;
    Eigen::Vector3d grad = k::cubic_spline_grad(xi, xj, h);
    CHECK(grad.norm() == Catch::Approx(0.0).margin(1e-18));
}

TEST_CASE("cubic_spline_grad: compact support", "[kernel][grad]")
{
    const double h = 0.5;
    Eigen::Vector3d xi(0.0, 0.0, 0.0);
    Eigen::Vector3d xj(1.5, 0.0, 0.0); // q = 3 > 2
    Eigen::Vector3d grad = k::cubic_spline_grad(xi, xj, h);
    CHECK(grad.norm() == Catch::Approx(0.0).margin(1e-18));
}

TEST_CASE("cubic_spline: compact support random", "[kernel][val]")
{
    const double h = 0.5;
    const double r_max = 5.0; // häufig > 2h
    std::mt19937 rng(12345); //seed for beeing reproducible
    std::uniform_real_distribution<double> dist(-r_max, r_max); //alternative to writing: "const double val = -r_max + (rng() / (double)rng.max()) * 2.0 * r_max"
    
    for (int i = 0; i < 100; ++i) {
        Eigen::Vector3d xi(dist(rng), dist(rng), dist(rng));
        Eigen::Vector3d xj(dist(rng), dist(rng), dist(rng));
        const double r = (xi-xj).norm();
        if (r > 2.0*h + 1e-12) { // outside support
            auto val = k::cubic_spline(xi, xj, h);
            CHECK(val == Catch::Approx(0.0).margin(1e-18));
        }
    }
}

TEST_CASE("cubic_spline_grad: compact support random", "[kernel][grad]")
{
    const double h = 0.5;
    const double r_max = 5.0; // häufig > 2h
    std::mt19937 rng(12345); //seed for beeing reproducible
    std::uniform_real_distribution<double> dist(-r_max, r_max); //alternative to writing: "const double val = -r_max + (rng() / (double)rng.max()) * 2.0 * r_max"
    
    for (int i = 0; i < 100; ++i) {
        Eigen::Vector3d xi(dist(rng), dist(rng), dist(rng));
        Eigen::Vector3d xj(dist(rng), dist(rng), dist(rng));
        const double r = (xi-xj).norm();
        if (r > 2.0*h + 1e-12) { // outside support
            auto grad = k::cubic_spline_grad(xi, xj, h);
            CHECK(grad.norm() == Catch::Approx(0.0).margin(1e-18));
        }
    }
}