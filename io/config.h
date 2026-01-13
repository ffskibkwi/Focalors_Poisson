#pragma once

#include <string>

class EnvironmentConfig
{
public:
    bool showCurrentStep = false; // Print current step
    bool showGmresRes    = false; // Print residual history of GMRES

    // Debug Output Control
    bool        debugMode      = false;
    std::string debugOutputDir = "./debug_output/";

    EnvironmentConfig() = default;
};

class TimeAdvancingConfig
{
public:
    double dt             = 0.0; // Time step
    double t_max          = 0.0; // Max time
    int    num_iterations = 0;   // Total number of time iterations
    int    corr_iter      = 1;   // Number of correction iterations per time step
    TimeAdvancingConfig() {};
    void set_dt(double in_dt);
    void set_t_max(double in_t_max);
    void set_num_iterations(int in_num_iterations);
    void set_corr_iter(int in_corr_iter);
};

class PhysicsConfig
{
public:
    double nu = 0.0;
    double Re = 0.0;

    // Non-Newtonian parameters
    // 0: Newtonian, 1: Power Law, 2: Carreau
    int model_type = 0;

    // Power Law
    double k = 0.0; // Consistency index

    // Carreau
    double mu_0   = 0.0; // Zero shear viscosity
    double mu_inf = 0.0; // Infinite shear viscosity
    double a      = 2.0; // Carreau parameter a
    double lambda = 0.0; // Relaxation time

    // Common
    double n      = 1.0;  // Power law index
    double mu_min = 0.0;  // Minimum viscosity limit
    double mu_max = 1e20; // Maximum viscosity limit

    PhysicsConfig() {}
    void set_nu(double in_nu);
    void set_Re(double in_Re);

    // Basic parameter setters
    void set_model_type(int in_model_type);
    void set_mu_min(double in_mu_min);
    void set_mu_max(double in_mu_max);

    // Power Law model setter (K and n)
    void set_power_law(double in_k, double in_n);
    // Power Law model setter (Re_PL and n)
    // Re_PL = rho * U^(2-n) * L^n / k
    // k = 1.0 / Re_PL (assuming rho=1, U=1, L=1)
    // Optional parameters in_mu_min and in_mu_max allow user to override default viscosity limits.
    // Use -1.0 as sentinel value to indicate "use default value".
    void set_power_law_dimensionless(double in_Re_PL, double in_n, double in_mu_min = -1.0, double in_mu_max = -1.0);

    // Carreau model setter (mu_0, mu_inf, a, lambda, n)
    void set_carreau(double in_mu_0, double in_mu_inf, double in_a, double in_lambda, double in_n);
    // Carreau model setter (Re_0, Re_inf, Wi, a, n)
    // Re_0 = rho * U * L / mu_0   => mu_0 = 1.0 / Re_0
    // Re_inf = rho * U * L / mu_inf => mu_inf = 1.0 / Re_inf
    // Wi = lambda * U / L        => lambda = Wi
    void set_carreau_dimensionless(double in_Re_0, double in_Re_inf, double in_Wi, double in_a, double in_n);
};
