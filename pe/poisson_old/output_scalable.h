#pragma once

class OutputScalable
{
public:
    void set_output_size(int nx, int ny)
    {
        output_nx = nx;
        output_ny = ny;
    }

    void set_output_size(int nx, int ny, int nz)
    {
        output_nx = nx;
        output_ny = ny;
        output_nz = nz;
    }

protected:
    int output_nx, output_ny, output_nz;
};