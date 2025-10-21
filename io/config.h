#pragma once

class EnvironmentConfig
{
public:
    bool showCurrentStep = false; // Print current step
    bool showGmresRes = false; // Print residual history of GMRES

    EnvironmentConfig() {};
};

class TimeAdvancingConfig
{
public:
    double dt = 0.0; //Time step
    double t_max = 0.0; //Max time
    int num_iterations = 0; //Total number of time iterations

    TimeAdvancingConfig() {};
    void set_dt(double in_dt);
    void set_t_max(double in_t_max);
    void set_num_iterations(int in_num_iterations);
};

class PhysicsConfig
{
public:
    double nu = 0.0;
    double Re = 0.0;

    PhysicsConfig() {}
    void set_nu(double in_nu);
    void set_Re(double in_Re);
};