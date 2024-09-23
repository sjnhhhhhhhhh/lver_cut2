#include <Open3D/Open3D.h>
#include <iostream>

int main() {
    // 创建简单的点云
    auto o3d_pcd = std::make_shared<open3d::geometry::PointCloud>();
    o3d_pcd->points_ = {
        {0.0, 0.0, 0.0},
        {1.0, 0.0, 0.0},
        {0.5, 1.0, 0.0},
        {0.5, 0.5, 1.0}
    };

    try {
        // 估计法向量
        o3d_pcd->EstimateNormals(open3d::geometry::KDTreeSearchParamHybrid(1.0, 30));
        o3d_pcd->OrientNormalsConsistentTangentPlane(100);

        // 执行泊松重建
        auto [o3d_mesh, densities] = open3d::geometry::TriangleMesh::CreateFromPointCloudPoisson(*o3d_pcd, 6);

        if (o3d_mesh->vertices_.empty() || o3d_mesh->triangles_.empty()) {
            std::cerr << "Error: Poisson reconstruction failed." << std::endl;
            return -1;
        }

        std::cout << "Poisson reconstruction succeeded." << std::endl;

    } catch (const std::exception& e) {
        std::cerr << "Exception during Poisson reconstruction: " << e.what() << std::endl;
        return -1;
    }

    return 0;
}
