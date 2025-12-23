#include "config.h"

void TimeAdvancingConfig::set_dt(double in_dt)
{
    dt = in_dt;
    if (t_max != 0.0 && num_iterations == 0)
    {
        num_iterations = static_cast<int>(t_max / dt);
    }
    else if (t_max == 0.0 && num_iterations != 0)
    {
        t_max = static_cast<double>(num_iterations) * dt;
    }
    else if (t_max != 0.0 && num_iterations != 0)
    {
        // We think total number of iterations is more important
        t_max = static_cast<double>(num_iterations) * dt;
    }
}

void TimeAdvancingConfig::set_t_max(double in_t_max)
{
    t_max = in_t_max;
    if (dt != 0.0)
    {
        num_iterations = static_cast<int>(t_max / dt);
    }
}

void TimeAdvancingConfig::set_num_iterations(int in_num_iterations)
{
    num_iterations = in_num_iterations;
    if (dt != 0.0)
    {
        t_max = static_cast<double>(num_iterations) * dt;
    }
}

void TimeAdvancingConfig::set_corr_iter(int in_corr_iter) { corr_iter = in_corr_iter; }

void PhysicsConfig::set_nu(double in_nu)
{
    nu = in_nu;
    Re = 1.0 / nu;
}

void PhysicsConfig::set_Re(double in_Re)
{
    Re = in_Re;
    nu = 1.0 / Re;
}

// Basic parameter setters
void PhysicsConfig::set_model_type(int in_model_type) { model_type = in_model_type; }

void PhysicsConfig::set_mu_min(double in_mu_min) { mu_min = in_mu_min; }

void PhysicsConfig::set_mu_max(double in_mu_max) { mu_max = in_mu_max; }

// Power Law model setter (K and n)
void PhysicsConfig::set_power_law(double in_k, double in_n)
{
    k = in_k;
    n = in_n;
}

void PhysicsConfig::set_power_law_dimensionless(double in_Re_PL, double in_n, double in_mu_min, double in_mu_max)
{
    // Re_PL = rho * U^(2-n) * L^n / k
    // k = 1.0 / Re_PL
    k = 1.0 / in_Re_PL;
    n = in_n;

    // For shear-thinning fluids (n <= 1), set viscosity limits to avoid numerical instability.
    // Power Law viscosity: mu = K * gamma_dot^(n-1)
    // When n < 1 (n-1 < 0):
    //   - As gamma_dot -> 0, mu -> infinity (causes numerical blow-up)
    //   - As gamma_dot -> infinity, mu -> 0 (may cause division issues)
    // Therefore, we need to clamp viscosity within reasonable bounds.
    if (in_n <= 1.0)
    {
        // mu_max: Upper limit to prevent viscosity from becoming too large at low shear rates
        // If user specified a valid value (>= 0), use it; otherwise use default (100 * k)
        if (in_mu_max >= 0.0)
        {
            mu_max = in_mu_max;
        }
        else
        {
            // Default: 100 times the reference viscosity (k = 1/Re_PL)
            mu_max = 10000.0 / in_Re_PL;
        }

        // mu_min: Lower limit to prevent viscosity from becoming too small at high shear rates
        // If user specified a valid value (>= 0), use it; otherwise use default (1e-4)
        if (in_mu_min >= 0.0)
        {
            mu_min = in_mu_min;
        }
        else
        {
            // Default: 0.01 times the reference viscosity (k = 1/Re_PL)
            mu_min = 0.01 / in_Re_PL;
        }
    }
}

// Carreau model setter (mu_0, mu_inf, a, lambda, n)
void PhysicsConfig::set_carreau(double in_mu_0, double in_mu_inf, double in_a, double in_lambda, double in_n)
{
    mu_0   = in_mu_0;
    mu_inf = in_mu_inf;
    a      = in_a;
    lambda = in_lambda;
    n      = in_n;
}

void PhysicsConfig::set_carreau_dimensionless(double in_Re_0, double in_Re_inf, double in_Wi, double in_a, double in_n)
{
    // Re_0 = rho * U * L / mu_0   => mu_0 = 1.0 / Re_0
    mu_0 = 1.0 / in_Re_0;

    // Re_inf = rho * U * L / mu_inf => mu_inf = 1.0 / Re_inf
    mu_inf = 1.0 / in_Re_inf;

    // Wi = lambda * U / L        => lambda = Wi
    lambda = in_Wi;

    a = in_a;
    n = in_n;
}
