#include "config.h"

void TimeAdvancingConfig::set_dt(double in_dt)
{
    dt = in_dt;
    if (t_max != 0.0 && num_iterations == 0)
    {
        num_iterations = static_cast<int>(t_max / dt);
    }else if(t_max == 0.0 && num_iterations != 0)
    {
        t_max = static_cast<double>(num_iterations) * dt;
    }else if(t_max != 0.0 && num_iterations != 0)
    {
        //We think total number of iterations is more important
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