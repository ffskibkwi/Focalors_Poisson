#include "base/domain/domain3d.h"
#include "base/domain/geometry3d.h"
#include "base/domain/geometry_tree.hpp"
#include "base/domain/variable3d.h"
#include "base/field/field3.h"

#include "base/location_boundary.h"

#include "pe/concat/concat_solver3d.h"

#include "io/config.h"
#include "io/csv_writer_3d.h"

int main(int argc, char* argv[])
{
    Geometry3D         geo_tee;
    EnvironmentConfig* env_config = new EnvironmentConfig();
    env_config->showGmresRes      = true;
    env_config->showCurrentStep   = true;

    // Grid parameters - 3D T型域
    // T2: 中心域 (x-y平面中心，z方向延伸)
    Domain3DUniform T2(10, 10, 10, 1.0, 1.0, 1.0, "T2");

    // T1: 左侧域 (x方向延伸)
    Domain3DUniform T1("T1");
    T1.set_nx(10);
    T1.set_lx(1.0);

    // T3: 右侧域 (x方向延伸)
    Domain3DUniform T3("T3");
    T3.set_nx(10);
    T3.set_lx(1.0);

    // T4: 下侧域 (y方向延伸)
    Domain3DUniform T4("T4");
    T4.set_nz(10);
    T4.set_lz(1.0);

    // Construct geometry
    geo_tee.add_domain(&T1);
    geo_tee.add_domain(&T2);
    geo_tee.add_domain(&T3);
    geo_tee.add_domain(&T4);

    // 3D连接：T型结构
    // T2为中心，左边连接T1，右边连接T3，下边连接T4
    geo_tee.connect(&T2, LocationType::Left, &T1);
    geo_tee.connect(&T2, LocationType::Right, &T3);
    geo_tee.connect(&T2, LocationType::Down, &T4);

    Variable3D v("v");
    v.set_geometry(geo_tee);

    // 创建field3并设置到Variable3D
    field3 v_T1("v_T1"), v_T2("v_T2"), v_T3("v_T3"), v_T4("v_T4");

    // 使用Variable3D的set_center_field方法
    v.set_center_field(&T1, v_T1);
    v.set_center_field(&T2, v_T2);
    v.set_center_field(&T3, v_T3);
    v.set_center_field(&T4, v_T4);

    // 设置边界条件
    // T2: 中心域，上边为Neumann边界
    v.set_boundary_type(&T2, LocationType::Front, PDEBoundaryType::Neumann);
    v.set_boundary_type(&T2, LocationType::Back, PDEBoundaryType::Neumann);
    v.set_boundary_type(&T2, LocationType::Up, PDEBoundaryType::Neumann);

    // T1: 左侧域，左、上、下、前、后为Neumann边界
    v.set_boundary_type(&T1, LocationType::Left, PDEBoundaryType::Neumann);
    v.set_boundary_type(&T1, LocationType::Front, PDEBoundaryType::Neumann);
    v.set_boundary_type(&T1, LocationType::Back, PDEBoundaryType::Neumann);
    v.set_boundary_type(&T1, LocationType::Up, PDEBoundaryType::Neumann);
    v.set_boundary_type(&T1, LocationType::Down, PDEBoundaryType::Neumann);

    // T3: 右侧域，右、上、下、前、后为Neumann边界
    v.set_boundary_type(&T3, LocationType::Right, PDEBoundaryType::Neumann);
    v.set_boundary_type(&T3, LocationType::Front, PDEBoundaryType::Neumann);
    v.set_boundary_type(&T3, LocationType::Back, PDEBoundaryType::Neumann);
    v.set_boundary_type(&T3, LocationType::Up, PDEBoundaryType::Neumann);
    v.set_boundary_type(&T3, LocationType::Down, PDEBoundaryType::Neumann);

    // T4: 下侧域，左、右、下、前、后为Neumann边界
    v.set_boundary_type(&T4, LocationType::Left, PDEBoundaryType::Neumann);
    v.set_boundary_type(&T4, LocationType::Right, PDEBoundaryType::Neumann);
    v.set_boundary_type(&T4, LocationType::Front, PDEBoundaryType::Neumann);
    v.set_boundary_type(&T4, LocationType::Back, PDEBoundaryType::Neumann);
    v.set_boundary_type(&T4, LocationType::Down, PDEBoundaryType::Neumann);

    // 设置初始值（在T4域的中心区域）
    for (int i = 3; i < v_T4.get_nx() - 3; i++)
    {
        for (int j = 3; j < v_T4.get_ny() - 3; j++)
        {
            for (int k = 3; k < v_T4.get_nz() - 3; k++)
            {
                v_T4(i, j, k) = 1.0;
            }
        }
    }

    ConcatPoissonSolver3D solver(&v, env_config);
    solver.solve();

    IO::write_csv(v_T1, "result/concat_new_3d/v_T1");
    IO::write_csv(v_T2, "result/concat_new_3d/v_T2");
    IO::write_csv(v_T3, "result/concat_new_3d/v_T3");
    IO::write_csv(v_T4, "result/concat_new_3d/v_T4");

    return 0;
}
