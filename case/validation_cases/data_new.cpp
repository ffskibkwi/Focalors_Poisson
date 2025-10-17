#include "base/domain/geometry2d.h"
#include "base/domain/domain2d.h"
#include "base/domain/variable.h"
#include "base/field/field2.h"
#include "base/domain/geometry_tree.hpp"

#include "base/location_boundary.h"

// #include "pe/concat/concat_solver2d.h" // 未在此测试中使用，避免引入外部依赖

inline const char* locationTypeToString(LocationType t)
    {
        switch (t)
        {
        case LocationType::Left:  return "Left";
        case LocationType::Right: return "Right";
        case LocationType::Down:  return "Down";
        case LocationType::Up:    return "Up";
        case LocationType::Front: return "Front";
        case LocationType::Back:  return "Back";
        default: return "Unknown";
        }
    }

int main(int argc, char* argv[])
{
    // 复杂几何测试：以 A 为中心的“星形”结构，测试单连通性检查与最优树生成
    //     C
    //     |
    // B -- A -- D
    //     |
    //     E

    // A 完整规格
    int    nx_A = 16;
    int    ny_A = 12;
    double lx_A = 1.6;
    double ly_A = 1.2;

    // 对于左右拼接(B、D)，预先设置 x 方向规格；y 方向由连接时对齐 A
    int    nx_B = 12; double lx_B = 1.2; // Left of A
    int    nx_D = 20; double lx_D = 2.0; // Right of A

    // 对于上下拼接(C、E)，预先设置 y 方向规格；x 方向由连接时对齐 A
    int    ny_C = 14; double ly_C = 1.4; // Up of A
    int    ny_E = 10; double ly_E = 1.0; // Down of A

    Geometry2D geo;
    Domain2DUniform A(nx_A, ny_A, lx_A, ly_A, "A");
    Domain2DUniform B("B");
    Domain2DUniform C("C");
    Domain2DUniform D("D");
    Domain2DUniform E("E");

    // 先设置将被连接方向的“垂直”规格，避免 profile 校验失败
    // 左/右相连：先给 B、D 设置 x 方向尺寸与节点
    B.set_nx(nx_B); B.set_lx(lx_B);
    D.set_nx(nx_D); D.set_lx(lx_D);
    // 上/下相连：先给 C、E 设置 y 方向尺寸与节点
    C.set_ny(ny_C); C.set_ly(ly_C);
    E.set_ny(ny_E); E.set_ly(ly_E);

    // 将域加入几何，并建立拼接关系
    geo.add_domain(A);
    geo.add_domain(B);
    geo.add_domain(C);
    geo.add_domain(D);
    geo.add_domain(E);

    geo.connect(A, LocationType::Left,  B);
    geo.connect(A, LocationType::Up,    C);
    geo.connect(A, LocationType::Right, D);
    geo.connect(A, LocationType::Down,  E);

    // 设置各域其余边界条件（连接的那一侧在 connect 中已设为 Dirichlet）
    A.set_boundary(LocationType::Left,  PDEBoundaryType::Dirichlet);
    A.set_boundary(LocationType::Right, PDEBoundaryType::Dirichlet);
    A.set_boundary(LocationType::Up,    PDEBoundaryType::Dirichlet);
    A.set_boundary(LocationType::Down,  PDEBoundaryType::Dirichlet);

    B.set_boundary(LocationType::Left,  PDEBoundaryType::Dirichlet);
    B.set_boundary(LocationType::Up,    PDEBoundaryType::Dirichlet);
    B.set_boundary(LocationType::Down,  PDEBoundaryType::Dirichlet);

    C.set_boundary(LocationType::Up,    PDEBoundaryType::Dirichlet);
    C.set_boundary(LocationType::Left,  PDEBoundaryType::Dirichlet);
    C.set_boundary(LocationType::Right, PDEBoundaryType::Dirichlet);

    D.set_boundary(LocationType::Right, PDEBoundaryType::Dirichlet);
    D.set_boundary(LocationType::Up,    PDEBoundaryType::Dirichlet);
    D.set_boundary(LocationType::Down,  PDEBoundaryType::Dirichlet);

    E.set_boundary(LocationType::Left,  PDEBoundaryType::Dirichlet);
    E.set_boundary(LocationType::Right, PDEBoundaryType::Dirichlet);
    E.set_boundary(LocationType::Down,  PDEBoundaryType::Dirichlet);

    // 1) 单连通性与主域检查
    try {
        geo.check();
        std::cout << "[OK] geometry check passed." << std::endl;
    } catch (const std::exception& ex) {
        std::cout << "[FAIL] geometry check threw: " << ex.what() << std::endl;
    }

    // 2) 构建最优树并打印（需要 geometry_tree.hpp 的 printTree）
    try {
        geo.solve_prepare();
        std::cout << "Geometry optimal tree:" << std::endl;
        TreeUtils::printTreeMap(geo.tree_root, geo.tree_map);
        if (geo.tree_root)
            std::cout << "tree root: " << geo.tree_root->name << std::endl;
    } catch (const std::exception& ex) {
        std::cout << "[FAIL] solve_prepare/printTree threw: " << ex.what() << std::endl;
    }

    // 3) 变量与场的绑定，确认映射
    Variable p("p");
    field2 p_A("p_A"), p_B("p_B"), p_C("p_C"), p_D("p_D"), p_E("p_E");
    p.set_geometry(geo);
    p.set_center_field(A, p_A);
    p.set_center_field(B, p_B);
    p.set_center_field(C, p_C);
    p.set_center_field(D, p_D);
    p.set_center_field(E, p_E);
    std::cout << "tree root of geometry: " << (geo.tree_root ? geo.tree_root->name : std::string("(null)")) << std::endl;
    std::cout << "fields bound: "
              << p_A.get_name() << ", " << p_B.get_name() << ", "
              << p_C.get_name() << ", " << p_D.get_name() << ", "
              << p_E.get_name() << std::endl;

    // 4) 负例一：非单连通（加入未连接的孤立子域）
    {
        Geometry2D geo_disconnected;
        Domain2DUniform X(8, 8, 0.8, 0.8, "X");
        Domain2DUniform Y("Y"); Y.set_nx(6); Y.set_lx(0.6); // 将与 X 相连
        Domain2DUniform Z(10, 10, 1.0, 1.0, "Z"); // 与任何域不连接，制造非单连通

        geo_disconnected.add_domain(X);
        geo_disconnected.add_domain(Y);
        geo_disconnected.add_domain(Z);
        geo_disconnected.connect(X, LocationType::Right, Y);

        X.set_boundary(LocationType::Left,  PDEBoundaryType::Dirichlet);
        X.set_boundary(LocationType::Up,    PDEBoundaryType::Dirichlet);
        X.set_boundary(LocationType::Down,  PDEBoundaryType::Dirichlet);
        Y.set_boundary(LocationType::Left,  PDEBoundaryType::Dirichlet);
        Y.set_boundary(LocationType::Up,    PDEBoundaryType::Dirichlet);
        Y.set_boundary(LocationType::Down,  PDEBoundaryType::Dirichlet);
        Z.set_boundary(LocationType::Left,  PDEBoundaryType::Dirichlet);
        Z.set_boundary(LocationType::Right, PDEBoundaryType::Dirichlet);
        Z.set_boundary(LocationType::Up,    PDEBoundaryType::Dirichlet);
        Z.set_boundary(LocationType::Down,  PDEBoundaryType::Dirichlet);

        try {
            geo_disconnected.check();
            std::cout << "[UNEXPECTED] disconnected geometry passed check." << std::endl;
        } catch (const std::exception& ex) {
            std::cout << "[EXPECTED FAIL] disconnected geometry check: " << ex.what() << std::endl;
        }
    }

    // 负例二：多主域（存在多个度>1 的节点）
    {
        Geometry2D geo_multi_main;
        Domain2DUniform P(6, 6, 0.6, 0.6, "P");
        Domain2DUniform Q("Q"); Q.set_nx(6); Q.set_lx(0.6);
        Domain2DUniform R("R"); R.set_nx(6); R.set_lx(0.6);
        Domain2DUniform S("S"); S.set_nx(6); S.set_lx(0.6);

        geo_multi_main.add_domain(P);
        geo_multi_main.add_domain(Q);
        geo_multi_main.add_domain(R);
        geo_multi_main.add_domain(S);

        // 构造 P - Q - R - S 的链路，Q 与 R 的度均为 2 -> 触发“多主域”
        geo_multi_main.connect(P, LocationType::Right, Q);
        geo_multi_main.connect(Q, LocationType::Right, R);
        geo_multi_main.connect(R, LocationType::Right, S);

        P.set_boundary(LocationType::Left,  PDEBoundaryType::Dirichlet);
        P.set_boundary(LocationType::Up,    PDEBoundaryType::Dirichlet);
        P.set_boundary(LocationType::Down,  PDEBoundaryType::Dirichlet);
        Q.set_boundary(LocationType::Up,    PDEBoundaryType::Dirichlet);
        Q.set_boundary(LocationType::Down,  PDEBoundaryType::Dirichlet);
        R.set_boundary(LocationType::Up,    PDEBoundaryType::Dirichlet);
        R.set_boundary(LocationType::Down,  PDEBoundaryType::Dirichlet);
        S.set_boundary(LocationType::Right, PDEBoundaryType::Dirichlet);
        S.set_boundary(LocationType::Up,    PDEBoundaryType::Dirichlet);
        S.set_boundary(LocationType::Down,  PDEBoundaryType::Dirichlet);

        try {
            geo_multi_main.check();
            std::cout << "[UNEXPECTED] multi-main geometry passed check." << std::endl;
        } catch (const std::exception& ex) {
            std::cout << "[EXPECTED FAIL] multi-main geometry check: " << ex.what() << std::endl;
        }
    }

    // 可选：后续联立求解
    // ConcatSolver2D solver(p);
    // solver.solve();

    // case1：检查可能有多 root 候选（路径 1-2-3-4-5-6，理论中心为 3 与 4）
    {
        std::cout << "case1: line 1-2-3-4-5-6" << std::endl;
        Geometry2D geo_line;
        Domain2DUniform L1(8, 8, 0.8, 0.8, "L1");
        Domain2DUniform L2("L2"); L2.set_nx(8); L2.set_lx(0.8);
        Domain2DUniform L3("L3"); L3.set_nx(8); L3.set_lx(0.8);
        Domain2DUniform L4("L4"); L4.set_nx(8); L4.set_lx(0.8);
        Domain2DUniform L5("L5"); L5.set_nx(8); L5.set_lx(0.8);
        Domain2DUniform L6("L6"); L6.set_nx(8); L6.set_lx(0.8);

        geo_line.add_domain(L1);
        geo_line.add_domain(L2);
        geo_line.add_domain(L3);
        geo_line.add_domain(L4);
        geo_line.add_domain(L5);
        geo_line.add_domain(L6);

        geo_line.connect(L1, LocationType::Right, L2);
        geo_line.connect(L2, LocationType::Right, L3);
        geo_line.connect(L3, LocationType::Right, L4);
        geo_line.connect(L4, LocationType::Right, L5);
        geo_line.connect(L5, LocationType::Right, L6);

        // 设置未连接侧的边界
        L1.set_boundary(LocationType::Left,  PDEBoundaryType::Dirichlet);
        L1.set_boundary(LocationType::Up,    PDEBoundaryType::Dirichlet);
        L1.set_boundary(LocationType::Down,  PDEBoundaryType::Dirichlet);

        L2.set_boundary(LocationType::Up,    PDEBoundaryType::Dirichlet);
        L2.set_boundary(LocationType::Down,  PDEBoundaryType::Dirichlet);
        L3.set_boundary(LocationType::Up,    PDEBoundaryType::Dirichlet);
        L3.set_boundary(LocationType::Down,  PDEBoundaryType::Dirichlet);
        L4.set_boundary(LocationType::Up,    PDEBoundaryType::Dirichlet);
        L4.set_boundary(LocationType::Down,  PDEBoundaryType::Dirichlet);
        L5.set_boundary(LocationType::Up,    PDEBoundaryType::Dirichlet);
        L5.set_boundary(LocationType::Down,  PDEBoundaryType::Dirichlet);

        L6.set_boundary(LocationType::Right, PDEBoundaryType::Dirichlet);
        L6.set_boundary(LocationType::Up,    PDEBoundaryType::Dirichlet);
        L6.set_boundary(LocationType::Down,  PDEBoundaryType::Dirichlet);

        try {
            geo_line.check();
            geo_line.solve_prepare();
            std::cout << "line optimal tree:" << std::endl;
            TreeUtils::printTreeMap(geo_line.tree_root, geo_line.tree_map);
            if (geo_line.tree_root)
                std::cout << "line tree root: " << geo_line.tree_root->name << std::endl;
        } catch (const std::exception& ex) {
            std::cout << "[FAIL] case1 threw: " << ex.what() << std::endl;
        }
    }

    // case2：复杂树状结构：横向 1-2-3；纵向 2-4-5；横向 5-6
    {
        std::cout << "case2: tee-like tree" << std::endl;
        Geometry2D geo_tee;
        Domain2DUniform T2(10, 10, 1.0, 1.0, "T2"); // 中心
        Domain2DUniform T1("T1"); T1.set_nx(10); T1.set_lx(1.0);
        Domain2DUniform T3("T3"); T3.set_nx(10); T3.set_lx(1.0);
        Domain2DUniform T4("T4"); T4.set_ny(10); T4.set_ly(1.0);
        Domain2DUniform T5("T5"); T5.set_ny(10); T5.set_ly(1.0);
        Domain2DUniform T6("T6"); T6.set_nx(10); T6.set_lx(1.0);

        geo_tee.add_domain(T1);
        geo_tee.add_domain(T2);
        geo_tee.add_domain(T3);
        geo_tee.add_domain(T4);
        geo_tee.add_domain(T5);
        geo_tee.add_domain(T6);

        // 横向 1-2-3
        geo_tee.connect(T2, LocationType::Left,  T1);
        geo_tee.connect(T2, LocationType::Right, T3);
        // 纵向 2-4-5
        geo_tee.connect(T2, LocationType::Down,  T4);
        geo_tee.connect(T4, LocationType::Down,  T5);
        // 横向 5-6
        geo_tee.connect(T5, LocationType::Right, T6);

        // 设置未连接侧的边界
        T2.set_boundary(LocationType::Up,    PDEBoundaryType::Dirichlet);

        T1.set_boundary(LocationType::Left,  PDEBoundaryType::Dirichlet);
        T1.set_boundary(LocationType::Up,    PDEBoundaryType::Dirichlet);
        T1.set_boundary(LocationType::Down,  PDEBoundaryType::Dirichlet);

        T3.set_boundary(LocationType::Right, PDEBoundaryType::Dirichlet);
        T3.set_boundary(LocationType::Up,    PDEBoundaryType::Dirichlet);
        T3.set_boundary(LocationType::Down,  PDEBoundaryType::Dirichlet);

        T4.set_boundary(LocationType::Left,  PDEBoundaryType::Dirichlet);
        T4.set_boundary(LocationType::Right, PDEBoundaryType::Dirichlet);

        T5.set_boundary(LocationType::Left,  PDEBoundaryType::Dirichlet);
        T5.set_boundary(LocationType::Down,  PDEBoundaryType::Dirichlet);

        T6.set_boundary(LocationType::Right, PDEBoundaryType::Dirichlet);
        T6.set_boundary(LocationType::Up,    PDEBoundaryType::Dirichlet);
        T6.set_boundary(LocationType::Down,  PDEBoundaryType::Dirichlet);

        try {
            geo_tee.check();
            geo_tee.solve_prepare();
            std::cout << "tee optimal tree:" << std::endl;
            TreeUtils::printTreeMap(geo_tee.tree_root, geo_tee.tree_map);
            
            if (geo_tee.tree_root)
                std::cout << "tee tree root: " << geo_tee.tree_root->name << std::endl;
            
            if (!geo_tee.parent_map.empty())
            {
                for (auto &[key, value] : geo_tee.parent_map)
                {
                    std::cout << "parent of " << key->name << " is " << value.second->name << " on " << locationTypeToString(value.first) << std::endl;
                }
                std::cout << std::endl;
            }
        } catch (const std::exception& ex) {
            std::cout << "[FAIL] case2 threw: " << ex.what() << std::endl;
        }
    }

    return 0;
}