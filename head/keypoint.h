#ifndef KEYPOINT_H
#define KEYPOINT_H

#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <Eigen/Dense>
#include <tuple>
#include <vector>

// 定义 PointWithLabel 结构体，用于标记点的归属
struct PointWithLabel {
    Eigen::Vector3d point;
    int label; // 标识点所属的组：1表示第一组，2表示第二组

    PointWithLabel(const Eigen::Vector3d& p, int l) : point(p), label(l) {}
};

// 提取 STL 文件中的点
std::vector<Eigen::Vector3d> extractpoints(vtkSmartPointer<vtkPolyData> data);

// 合并两个点集并为每个点打上标签
std::vector<PointWithLabel> mergePointSets(const std::vector<Eigen::Vector3d>& set1, const std::vector<Eigen::Vector3d>& set2);

// 使用体素栅格滤波器对点集进行下采样
void downsamplePoints(std::vector<Eigen::Vector3d>& points, double leafSize);

// 拟合二次多项式表面：z = ax^2 + by^2 + cxy + dx + ey + f
Eigen::VectorXd fitPolynomialSurface(const std::vector<Eigen::Vector3d>& points);

// 将 Eigen::Vector3d 点集转换为 VTK 格式的 PolyData
vtkSmartPointer<vtkPolyData> vec_to_poly(std::vector<Eigen::Vector3d> points);

// 可视化拟合的多项式表面
vtkSmartPointer<vtkPolyData> visualizePolynomialSurface(const Eigen::VectorXd& beta, double bounds[6], int resolution);

// 根据指定边界过滤点集
std::vector<Eigen::Vector3d> Filter_points(std::vector<Eigen::Vector3d> points, double bounds[6]);

// 计算两个边界框的最大范围
void max_bounds(double bounds1[6], double bounds2[6], double mbounds[6]);

std::vector<Eigen::Vector3d> uniformSampling(const std::vector<Eigen::Vector3d>& points, size_t num_samples);

// 从两个 VTK 点云中查找关键点，并返回过滤后的点集和 VTK 格式的点集
std::tuple<vtkSmartPointer<vtkPolyData>, std::vector<Eigen::Vector3d>> find_points(vtkSmartPointer<vtkPolyData> input1, vtkSmartPointer<vtkPolyData> input2);

// 函数：结合 vtkCleanPolyData 和 vtkVoxelGrid 进行点云去噪
vtkSmartPointer<vtkPolyData> denoisePointCloud(vtkSmartPointer<vtkPolyData> inputPolyData);



#endif // KEYPOINT_H
