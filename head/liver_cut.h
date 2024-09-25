#ifndef LIVER_CUT_H
#define LIVER_CUT_H
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
#include <vtkImageData.h>
#include <vtkMarchingCubes.h>
#include <vtkDoubleArray.h>
#include <vtkCellArray.h>
#include <vtkProperty.h>
#include <vtkSurfaceReconstructionFilter.h>
#include <vtkContourFilter.h>
#include <vtkImplicitPolyDataDistance.h>
#include <vtkClipPolyData.h>
#include <vtkSampleFunction.h>
#include <vtkAppendPolyData.h>
#include <vtkPCANormalEstimation.h>
#include <vtkSignedDistance.h>
#include <vtkSurfaceReconstructionFilter.h>
#include <vtkCleanPolyData.h>
#include <vtkVoxelGrid.h>

#include <iostream>
#include <vector>
#include <cmath>
#include <tuple>
#include <algorithm>
#include <stdexcept> 


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
#include <chrono>
#include <Eigen/Dense>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Alpha_shape_3.h>


#include <pcl/point_types.h>
#include <pcl/filters/voxel_grid.h>

typedef CGAL::Simple_cartesian<double> K;
typedef K::Point_3 Point_3;
typedef CGAL::Delaunay_triangulation_3<K> Delaunay;
typedef K::Segment_3 Segment_3;
typedef K::Plane_3 Plane_3;
typedef CGAL::Alpha_shape_3<K> Alpha_shape;


struct PointWithLabel {
    Eigen::Vector3d point;
    int label; // 标识点所属的组：1表示第一组，2表示第二组
};

// 将点集转换为 vtkImageData
vtkSmartPointer<vtkImageData> convertPointsToImageData(const std::vector<Eigen::Vector3d>& points, double voxelSize, int gridSize[3]);

// 提取STL文件中的点
std::vector<Eigen::Vector3d> extractpoints(vtkSmartPointer<vtkPolyData> data);

// 合并两个点集，并标记点的归属
std::vector<PointWithLabel> mergePointSets(const std::vector<Eigen::Vector3d>& set1, const std::vector<Eigen::Vector3d>& set2);

// 将 std::vector<Eigen::Vector3d> 转换为 vtkPolyData
vtkSmartPointer<vtkPolyData> vec_to_poly(std::vector<Eigen::Vector3d> points);

// 过滤掉不在范围内的点集
std::vector<Eigen::Vector3d> Filter_points(std::vector<Eigen::Vector3d> points, double bounds[6]);

// 计算两个边界的最大范围
void max_bounds(double bounds1[6], double bounds2[6], double mbounds[6]);

// 可变参数模板函数，用于合并任意数量的 vtkPolyData 对象
template <typename... PolyDataTypes>
vtkSmartPointer<vtkPolyData> mergePolyData(const PolyDataTypes&... inputs) {
    vtkSmartPointer<vtkAppendPolyData> appendFilter = vtkSmartPointer<vtkAppendPolyData>::New();

    // Lambda 表达式将每个输入添加到 appendFilter
    auto addInput = [&](const vtkSmartPointer<vtkPolyData>& input) {
        appendFilter->AddInputData(input);
    };

    // 使用折叠表达式遍历所有传入的 vtkPolyData 对象
    (addInput(inputs), ...);

    appendFilter->Update();

    vtkSmartPointer<vtkPolyData> mergedPolyData = vtkSmartPointer<vtkPolyData>::New();
    mergedPolyData->ShallowCopy(appendFilter->GetOutput());

    return mergedPolyData;
}

// 对点集进行均匀采样
std::vector<Eigen::Vector3d> uniformSampling(const std::vector<Eigen::Vector3d>& points, size_t num_samples);

// 计算点集的边界
void calculateBounds(const std::vector<Eigen::Vector3d>& points, double bounds[6]);

// 根据输入生成隐式平面
vtkSmartPointer<vtkImplicitPolyDataDistance> plane_generator(vtkSmartPointer<vtkPolyData> Poly_vtk_points);

// 使用隐式平面切割肝脏模型
vtkSmartPointer<vtkPolyData> cut_run(vtkSmartPointer<vtkImplicitPolyDataDistance> implicitFunction, vtkSmartPointer<vtkPolyData> liver);

// 结合 vtkCleanPolyData 和 vtkVoxelGrid 进行点云去噪
vtkSmartPointer<vtkPolyData> denoisePointCloud(vtkSmartPointer<vtkPolyData> inputPolyData, double cleanTolerance, double leafSizeX, double leafSizeY, double leafSizeZ);

// 使用隐式平面切割肝脏模型，返回切割后的主要部分和被切下的小块
std::tuple<vtkSmartPointer<vtkPolyData>, vtkSmartPointer<vtkPolyData>> cut_run2(vtkSmartPointer<vtkImplicitPolyDataDistance> implicitFunction, vtkSmartPointer<vtkPolyData> liver);

// 找到两个输入模型中的点并进行处理，返回去噪后的点和过滤点
std::tuple<vtkSmartPointer<vtkPolyData>, std::vector<Eigen::Vector3d>> find_points(vtkSmartPointer<vtkPolyData> input1, vtkSmartPointer<vtkPolyData> input2);

#endif // LIVER_CUT_H
