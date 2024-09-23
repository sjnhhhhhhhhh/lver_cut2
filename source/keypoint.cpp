#include <vtkSmartPointer.h>
#include <vtkRendererCollection.h>
#include <vtkPointPicker.h>
#include <vtkSphereSource.h>
#include <vtkPolyDataMapper.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkActor.h>
#include <vtkInteractorStyleTrackballCamera.h>
#include <vtkObjectFactory.h>
#include <vtkProperty.h>
#include <vtkSTLReader.h>
#include <vtkCylinderSource.h>
#include <vtkCubeAxesActor.h>
#include <vtkAxesActor.h>
#include <vtkOrientationMarkerWidget.h>
#include <vtkAutoInit.h>
#include <vtkPlaneSource.h>
#include <vtkDataSetMapper.h>
#include <vtkSTLWriter.h>
#include <vtkGeometryFilter.h>
#include <vtkTriangle.h>
#include <vtkVertexGlyphFilter.h>
#include <vtkPolyDataWriter.h>
#include <vtkPlane.h>
#include <vtkCutter.h>
#include <vtkLineSource.h>
#include <vtkConvexHull2D.h>
#include <vtkDelaunay3D.h>
#include <vtkGeometryFilter.h>
#include <vtkConvexPointSet.h>
#include <vtkUnstructuredGrid.h>
#include <vtkCommand.h>
#include <vtkCamera.h>
#include <vtkVoxelGrid.h>

#include <iostream>
#include <vector>
#include <cmath>
#include <tuple>
#include <chrono>
#include <unordered_map>

#include <Eigen/Dense>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Delaunay_triangulation_3.h>

#include <pcl/point_types.h>
#include <pcl/filters/voxel_grid.h>

#include "keypoint.h"

using namespace std;

typedef CGAL::Simple_cartesian<double> K;
typedef K::Point_3 Point_3;
typedef CGAL::Delaunay_triangulation_3<K> Delaunay;
typedef K::Segment_3 Segment_3;
typedef K::Plane_3 Plane_3;

/*struct PointWithLabel {
    Eigen::Vector3d point;
    int label; // 标识点所属的组：1表示第一组，2表示第二组
};*/

// 提取STL文件中的点
std::vector<Eigen::Vector3d> extractpoints(vtkSmartPointer<vtkPolyData> data) {
    vtkSmartPointer<vtkPoints> points = data->GetPoints();

    std::vector<Eigen::Vector3d> pointsvector;
    vtkIdType num = points->GetNumberOfPoints();
    for (vtkIdType i = 0; i < num; i++) {
        double temp[3];
        points->GetPoint(i, temp);
        Eigen::Vector3d p(temp[0], temp[1], temp[2]);
        pointsvector.push_back(p);
    }
    return pointsvector;
}

std::vector<PointWithLabel> mergePointSets(const std::vector<Eigen::Vector3d>& set1, const std::vector<Eigen::Vector3d>& set2) {
    std::vector<PointWithLabel> mergedPoints;

    for (const auto& point : set1) {
        mergedPoints.push_back(PointWithLabel(point, 1));  // 使用显式构造函数
    }
    for (const auto& point : set2) {
        mergedPoints.push_back(PointWithLabel(point, 2));  // 使用显式构造函数
    }

    return mergedPoints;
}

/*
std::vector<PointWithLabel> mergePointSets(const std::vector<Eigen::Vector3d>& set1, const std::vector<Eigen::Vector3d>& set2) {
    std::vector<PointWithLabel> mergedPoints;

    for (const auto& point : set1) {
        mergedPoints.push_back({point, 1});
    }
    for (const auto& point : set2) {
        mergedPoints.push_back({point, 2});
    }

    return mergedPoints;
}*/

void downsamplePoints(std::vector<Eigen::Vector3d>& points, double leafSize) {
    // 将 Eigen::Vector3d 点集转换为 PCL 的点云格式
    pcl::PointCloud<pcl::PointXYZ>::Ptr cloud(new pcl::PointCloud<pcl::PointXYZ>);
    for (const auto& p : points) {
        cloud->points.push_back(pcl::PointXYZ(p[0], p[1], p[2]));
    }
    
    // 创建体素栅格滤波器
    pcl::VoxelGrid<pcl::PointXYZ> sor;
    sor.setInputCloud(cloud);
    sor.setLeafSize(leafSize, leafSize, leafSize);

    // 执行滤波
    pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_filtered(new pcl::PointCloud<pcl::PointXYZ>);
    sor.filter(*cloud_filtered);

    // 将滤波后的点云转换回 Eigen::Vector3d
    points.clear();
    for (const auto& p : cloud_filtered->points) {
        points.push_back(Eigen::Vector3d(p.x, p.y, p.z));
    }
}

// 二次多项式拟合：拟合二次曲面 z = ax^2 + by^2 + cxy + dx + ey + f
Eigen::VectorXd fitPolynomialSurface(const std::vector<Eigen::Vector3d>& points) {
    int num_points = points.size();

    // 构建设计矩阵 X 和响应向量 Z
    Eigen::MatrixXd X(num_points, 6);  // 6 个多项式系数: a, b, c, d, e, f
    Eigen::VectorXd Z(num_points);     // z 值

    for (int i = 0; i < num_points; ++i) {
        double x = points[i][0];
        double y = points[i][1];
        double z = points[i][2];

        // 多项式系数的设计矩阵 X
        X(i, 0) = x * x;  // x^2
        X(i, 1) = y * y;  // y^2
        X(i, 2) = x * y;  // xy
        X(i, 3) = x;      // x
        X(i, 4) = y;      // y
        X(i, 5) = 1.0;    // 常数项

        // z 值
        Z(i) = z;
    }

    // 求解最小二乘问题 X * beta = Z，得到拟合系数 beta
    Eigen::VectorXd beta = X.colPivHouseholderQr().solve(Z);

    return beta;  // 返回拟合的系数: [a, b, c, d, e, f]
}

vtkSmartPointer<vtkPolyData> vec_to_poly(std::vector<Eigen::Vector3d> points){
    // 创建 VTK 点集
    vtkSmartPointer<vtkPoints> vtk_points = vtkSmartPointer<vtkPoints>::New();
    for (const auto& point : points) {
        vtk_points->InsertNextPoint(point[0], point[1], point[2]);
    }

    // 创建 PolyData 并设置点集
    vtkSmartPointer<vtkPolyData> polyData = vtkSmartPointer<vtkPolyData>::New();
    polyData->SetPoints(vtk_points);

    return polyData;
}

// 可视化拟合的多项式曲面
vtkSmartPointer<vtkPolyData> visualizePolynomialSurface(const Eigen::VectorXd& beta, double bounds[6], int resolution) {
    // 1. 创建 VTK 点集
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkCellArray> triangles = vtkSmartPointer<vtkCellArray>::New();

    double xMin = bounds[0];
    double xMax = bounds[1];
    double yMin = bounds[2];
    double yMax = bounds[3];
    // 网格划分
    double dx = (xMax - xMin) / (resolution - 1);
    double dy = (yMax - yMin) / (resolution - 1);

    std::vector<std::vector<vtkIdType>> pointIds(resolution, std::vector<vtkIdType>(resolution));

    // 2. 生成网格点并计算 z 值
    for (int i = 0; i < resolution; ++i) {
        for (int j = 0; j < resolution; ++j) {
            double x = xMin + i * dx;
            double y = yMin + j * dy;

            // 根据拟合的多项式系数计算 z 值
            double z = beta[0] * x * x + beta[1] * y * y + beta[2] * x * y + beta[3] * x + beta[4] * y + beta[5];

            // 插入点到 VTK 点集中
            pointIds[i][j] = points->InsertNextPoint(x, y, z);
        }
    }

    // 3. 生成三角形网格
    for (int i = 0; i < resolution - 1; ++i) {
        for (int j = 0; j < resolution - 1; ++j) {
            vtkSmartPointer<vtkTriangle> triangle1 = vtkSmartPointer<vtkTriangle>::New();
            triangle1->GetPointIds()->SetId(0, pointIds[i][j]);
            triangle1->GetPointIds()->SetId(1, pointIds[i + 1][j]);
            triangle1->GetPointIds()->SetId(2, pointIds[i][j + 1]);

            vtkSmartPointer<vtkTriangle> triangle2 = vtkSmartPointer<vtkTriangle>::New();
            triangle2->GetPointIds()->SetId(0, pointIds[i + 1][j + 1]);
            triangle2->GetPointIds()->SetId(1, pointIds[i][j + 1]);
            triangle2->GetPointIds()->SetId(2, pointIds[i + 1][j]);

            triangles->InsertNextCell(triangle1);
            triangles->InsertNextCell(triangle2);
        }
    }

    // 4. 创建 PolyData 并设置点集和三角形网格
    vtkSmartPointer<vtkPolyData> polyData = vtkSmartPointer<vtkPolyData>::New();
    polyData->SetPoints(points);
    polyData->SetPolys(triangles);
    return polyData;
}

std::vector<Eigen::Vector3d> Filter_points(std::vector<Eigen::Vector3d> points,double bounds[6]){
    std::vector<Eigen::Vector3d> filteredPoints;

    for (const auto& point : points) {
        // 检查点是否在 x, y, z 的范围内
        if (point[0] >= bounds[0] && point[0] <= bounds[1] &&  // x 轴范围
            point[1] >= bounds[2] && point[1] <= bounds[3] &&  // y 轴范围
            point[2] >= bounds[4] && point[2] <= bounds[5]) {  // z 轴范围
            filteredPoints.push_back(point);
        }
    }

    return filteredPoints;  // 返回过滤后的点集
} 

void max_bounds(double bounds1[6], double bounds2[6], double mbounds[6]) {
    // 计算最小的x, y, z边界
    mbounds[0] = (bounds1[0] < bounds2[0]) ? bounds1[0] : bounds2[0];
    mbounds[2] = (bounds1[2] < bounds2[2]) ? bounds1[2] : bounds2[2];
    mbounds[4] = (bounds1[4] < bounds2[4]) ? bounds1[4] : bounds2[4];

    // 计算最大的x, y, z边界
    mbounds[1] = (bounds1[1] > bounds2[1]) ? bounds1[1] : bounds2[1];
    mbounds[3] = (bounds1[3] > bounds2[3]) ? bounds1[3] : bounds2[3];
    mbounds[5] = (bounds1[5] > bounds2[5]) ? bounds1[5] : bounds2[5];
}

std::vector<Eigen::Vector3d> uniformSampling(const std::vector<Eigen::Vector3d>& points, size_t num_samples){
    std::vector<Eigen::Vector3d> sampled_points;
    size_t step = points.size() / num_samples;

    for (size_t i = 0; i < points.size(); i += step) {
        sampled_points.push_back(points[i]);
    }

    return sampled_points;    
}

// 函数：结合 vtkCleanPolyData 和 vtkVoxelGrid 进行点云去噪
vtkSmartPointer<vtkPolyData> denoisePointCloud(vtkSmartPointer<vtkPolyData> inputPolyData) {
    vtkSmartPointer<vtkVoxelGrid> voxelGrid = vtkSmartPointer<vtkVoxelGrid>::New();
    voxelGrid->SetInputData(inputPolyData);
    voxelGrid->SetLeafSize(1, 1, 1); // 设置体素大小
    voxelGrid->Update();
    
    vtkSmartPointer<vtkPolyData> denoisedPolyData = voxelGrid->GetOutput();
    return denoisedPolyData;
}

std::tuple<vtkSmartPointer<vtkPolyData>,std::vector<Eigen::Vector3d>>  find_points(vtkSmartPointer<vtkPolyData> input1, vtkSmartPointer<vtkPolyData> input2){
    // 获取最大边界，
    double boundsv1[6];
    double boundsv2[6];
    input1->GetBounds(boundsv1);
    //cout<<"r1v1 range:\n"<<"x:("<<boundsv1[0]<<","<<boundsv1[1]<<")\n";
    //cout<<"y:("<<boundsv1[2]<<","<<boundsv1[3]<<")\n";
    input2->GetBounds(boundsv2);
    //cout<<"r1v1 range:\n"<<"x:("<<boundsv2[0]<<","<<boundsv2[1]<<")\n";
    //cout<<"y:("<<boundsv2[2]<<","<<boundsv2[3]<<")\n";
    double mbounds[6];
    max_bounds(boundsv1,boundsv2,mbounds);

    // 合并点集，准备德拉内三角化
    auto pointsvector1 = extractpoints(input1);
    auto pointsvector2 = extractpoints(input2);
    // 合并点集并标记点的归属
    auto points = mergePointSets(pointsvector1, pointsvector2);
    // 切割平面的质心和法向量
    Eigen::Vector3d centroid(0, 0, 0);
    Eigen::Vector3d normal;

    // 创建 CGAL Point_3 的向量并构建哈希表
    std::vector<Point_3> cgal_points;
    std::unordered_map<Point_3, int> point_label_map;
    for (const auto& pointWithLabel : points) {
        Point_3 p(pointWithLabel.point[0], pointWithLabel.point[1], pointWithLabel.point[2]);
        cgal_points.push_back(p);
        point_label_map[p] = pointWithLabel.label;  // 哈希表存储点与标签的对应关系
    }

    // 生成 Delaunay 三角化
    Delaunay dt;
    dt.insert(cgal_points.begin(), cgal_points.end());

    std::vector<Eigen::Vector3d> points_for_fitting;  // 用于拟合的点集

    // 遍历 Delaunay 四面体（Cell），提取位于不同点集之间的面
    for (auto cit = dt.finite_cells_begin(); cit != dt.finite_cells_end(); ++cit) {
        Delaunay::Cell_handle cell = cit;

        // 遍历该四面体的每个面（Facet）
        for (int i = 0; i < 4; ++i) {
            // 获取面上三个顶点的标签
            int label1 = point_label_map[cell->vertex((i + 1) % 4)->point()];
            int label2 = point_label_map[cell->vertex((i + 2) % 4)->point()];
            int label3 = point_label_map[cell->vertex((i + 3) % 4)->point()];

            // 如果这三个顶点中的至少两个点来自不同的点集，则认为这是一个分割面
            if ((label1 != label2) || (label1 != label3)) {
                // 构建 Delaunay 面（Facet），并提取 Voronoi 对偶面
                Delaunay::Facet facet(cell, i);
                CGAL::Object voronoi_object = dt.dual(facet);

                // 尝试将 Voronoi 对偶对象转换为线段
                Segment_3 voronoi_segment;
                if (CGAL::assign(voronoi_segment, voronoi_object)) {
                    // 将 Voronoi 线段的两个端点加入拟合点集中
                    points_for_fitting.push_back(Eigen::Vector3d(voronoi_segment.source().x(), voronoi_segment.source().y(), voronoi_segment.source().z()));
                    points_for_fitting.push_back(Eigen::Vector3d(voronoi_segment.target().x(), voronoi_segment.target().y(), voronoi_segment.target().z()));
                }
            }
        }
    }
    //size = 32404
    //cout<<"size = "<<points_for_fitting.size()<<"\n";

    // 均匀稀疏化
    points_for_fitting = uniformSampling(points_for_fitting,10000);
    cout<<"size_processed = "<<points_for_fitting.size()<<"\n";
    
    auto filtered_points = Filter_points(points_for_fitting,mbounds);
    cout<<"size_processed = "<<filtered_points.size()<<"\n";

    vtkSmartPointer<vtkPoints> vtk_points = vtkSmartPointer<vtkPoints>::New();
    for (const auto& point : filtered_points) {
        vtk_points->InsertNextPoint(point[0], point[1], point[2]);
    }

    // 2. 创建 PolyData 并设置点集
    vtkSmartPointer<vtkPolyData> polyData_points = vtkSmartPointer<vtkPolyData>::New();
    polyData_points->SetPoints(vtk_points);

    auto vtk_points_denoise = denoisePointCloud(polyData_points);

    return std::make_tuple(vtk_points_denoise,filtered_points);
}