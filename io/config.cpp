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

// Carreau model setter (mu_0, mu_inf, a, lambda, n)
void PhysicsConfig::set_carreau(double in_mu_0, double in_mu_inf, double in_a, double in_lambda, double in_n)
{
    mu_0   = in_mu_0;
    mu_inf = in_mu_inf;
    a      = in_a;
    lambda = in_lambda;
    n      = in_n;
}