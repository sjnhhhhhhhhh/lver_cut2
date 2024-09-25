#include "liver_cut.h"
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
using namespace std;
// 将你的 filtered_points 转换为 vtkImageData
vtkSmartPointer<vtkImageData> convertPointsToImageData(const std::vector<Eigen::Vector3d>& points, double voxelSize, int gridSize[3]) {
    vtkSmartPointer<vtkImageData> imageData = vtkSmartPointer<vtkImageData>::New();

    // 设置体数据大小，体素间隔，和范围
    imageData->SetDimensions(gridSize[0], gridSize[1], gridSize[2]);
    imageData->SetSpacing(voxelSize, voxelSize, voxelSize);
    imageData->AllocateScalars(VTK_DOUBLE, 1);

    // 初始化体数据为0
    for (int z = 0; z < gridSize[2]; ++z) {
        for (int y = 0; y < gridSize[1]; ++y) {
            for (int x = 0; x < gridSize[0]; ++x) {
                double* pixel = static_cast<double*>(imageData->GetScalarPointer(x, y, z));
                *pixel = 0.0;
            }
        }
    }

    // 将 filtered_points 中的点设置为1，表示这些体素为表面
    for (const auto& point : points) {
        int x = static_cast<int>(point[0] / voxelSize);
        int y = static_cast<int>(point[1] / voxelSize);
        int z = static_cast<int>(point[2] / voxelSize);

        if (x >= 0 && x < gridSize[0] && y >= 0 && y < gridSize[1] && z >= 0 && z < gridSize[2]) {
            double* pixel = static_cast<double*>(imageData->GetScalarPointer(x, y, z));
            *pixel = 1.0; // 标记体素为非零值，代表物体表面
        }
    }

    return imageData;
}

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
        mergedPoints.push_back({point, 1});
    }
    for (const auto& point : set2) {
        mergedPoints.push_back({point, 2});
    }

    return mergedPoints;
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

// 可变参数模板函数，用于合并任意数量的 vtkPolyData 对象
/*template <typename... PolyDataTypes>
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
}*/

std::vector<Eigen::Vector3d> uniformSampling(const std::vector<Eigen::Vector3d>& points, size_t num_samples) {
    std::vector<Eigen::Vector3d> sampled_points;
    size_t step = points.size() / num_samples;

    for (size_t i = 0; i < points.size(); i += step) {
        sampled_points.push_back(points[i]);
    }

    return sampled_points;
}

void calculateBounds(const std::vector<Eigen::Vector3d>& points, double bounds[6]) {
    if (points.empty()) return;

    bounds[0] = bounds[1] = points[0][0];
    bounds[2] = bounds[3] = points[0][1];
    bounds[4] = bounds[5] = points[0][2];

    for (const auto& point : points) {
        if (point[0] < bounds[0]) bounds[0] = point[0];
        if (point[0] > bounds[1]) bounds[1] = point[0];
        if (point[1] < bounds[2]) bounds[2] = point[1];
        if (point[1] > bounds[3]) bounds[3] = point[1];
        if (point[2] < bounds[4]) bounds[4] = point[2];
        if (point[2] > bounds[5]) bounds[5] = point[2];
    }
}

vtkSmartPointer<vtkImplicitPolyDataDistance> plane_generator(vtkSmartPointer<vtkPolyData> Poly_vtk_points){


    vtkSmartPointer<vtkSurfaceReconstructionFilter> surfaceReconstruction = vtkSmartPointer<vtkSurfaceReconstructionFilter>::New();
    surfaceReconstruction->SetInputData(Poly_vtk_points);
    surfaceReconstruction->Update();

    vtkSmartPointer<vtkContourFilter> contourFilter = vtkSmartPointer<vtkContourFilter>::New();
    contourFilter->SetInputConnection(surfaceReconstruction->GetOutputPort());
    contourFilter->SetValue(0, 0.0);
    contourFilter->Update();

    // 获取拟合曲面的 vtkPolyData
    vtkSmartPointer<vtkPolyData> fittedSurface = contourFilter->GetOutput();

    // 创建一个隐式函数
    vtkSmartPointer<vtkImplicitPolyDataDistance> implicitFunction = vtkSmartPointer<vtkImplicitPolyDataDistance>::New();
    implicitFunction->SetInput(fittedSurface);


    return implicitFunction;
}

vtkSmartPointer<vtkPolyData> cut_run(vtkSmartPointer<vtkImplicitPolyDataDistance> implicitFunction, vtkSmartPointer<vtkPolyData> liver){

    vtkSmartPointer<vtkClipPolyData> clipper = vtkSmartPointer<vtkClipPolyData>::New();
    clipper->SetInputData(liver);
    clipper->SetClipFunction(implicitFunction);
    clipper->SetInsideOut(true); // 保留切下来的部分（小块）
    clipper->Update();

    vtkSmartPointer<vtkPolyData> clippedModel = clipper->GetOutput();
    return clippedModel;
}

// 函数：结合 vtkCleanPolyData 和 vtkVoxelGrid 进行点云去噪
vtkSmartPointer<vtkPolyData> denoisePointCloud(vtkSmartPointer<vtkPolyData> inputPolyData, double cleanTolerance, double leafSizeX, double leafSizeY, double leafSizeZ) {
    vtkSmartPointer<vtkVoxelGrid> voxelGrid = vtkSmartPointer<vtkVoxelGrid>::New();
    voxelGrid->SetInputData(inputPolyData);
    voxelGrid->SetLeafSize(leafSizeX, leafSizeY, leafSizeZ); // 设置体素大小
    voxelGrid->Update();
    
    vtkSmartPointer<vtkPolyData> denoisedPolyData = voxelGrid->GetOutput();
    return denoisedPolyData;
}

std::tuple<vtkSmartPointer<vtkPolyData>,vtkSmartPointer<vtkPolyData>>cut_run2(vtkSmartPointer<vtkImplicitPolyDataDistance> implicitFunction, vtkSmartPointer<vtkPolyData> liver){
    // 1. 获取裁剪后的主要部分
    vtkSmartPointer<vtkClipPolyData> clipperMain = vtkSmartPointer<vtkClipPolyData>::New();
    clipperMain->SetInputData(liver);
    clipperMain->SetClipFunction(implicitFunction);
    clipperMain->SetInsideOut(false); // 保留外部部分
    clipperMain->Update();
    vtkSmartPointer<vtkPolyData> mainModel = clipperMain->GetOutput();

    // 2. 获取被切下的小块
    vtkSmartPointer<vtkClipPolyData> clipperCut = vtkSmartPointer<vtkClipPolyData>::New();
    clipperCut->SetInputData(liver);
    clipperCut->SetClipFunction(implicitFunction);
    clipperCut->SetInsideOut(true); // 保留内部部分
    clipperCut->Update();
    vtkSmartPointer<vtkPolyData> cutOffPiece = clipperCut->GetOutput();

    return std::make_tuple(mainModel,cutOffPiece);

}

std::tuple<vtkSmartPointer<vtkPolyData>,std::vector<Eigen::Vector3d>>  find_points(vtkSmartPointer<vtkPolyData> input1, vtkSmartPointer<vtkPolyData> input2){
    
    
    double boundsv1[6];
    double boundsv2[6];

    input1->GetBounds(boundsv1);
    cout<<"r1v1 range:\n"<<"x:("<<boundsv1[0]<<","<<boundsv1[1]<<")\n";
    cout<<"y:("<<boundsv1[2]<<","<<boundsv1[3]<<")\n";
    cout<<"z:("<<boundsv1[4]<<","<<boundsv1[5]<<")\n";

    input2->GetBounds(boundsv2);
    cout<<"r1v1 range:\n"<<"x:("<<boundsv2[0]<<","<<boundsv2[1]<<")\n";
    cout<<"y:("<<boundsv2[2]<<","<<boundsv2[3]<<")\n";
    cout<<"z:("<<boundsv2[4]<<","<<boundsv2[5]<<")\n";

    double mbounds[6];
    max_bounds(boundsv1,boundsv2,mbounds);


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
    cout<<"size = "<<points_for_fitting.size()<<"\n";

    //downsamplePoints(points_for_fitting,200+i*10);
    //cout<<"size_processed = "<<points_for_fitting.size()<<"\n";

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
    double cleanTolerance = 1; // 合并容差
    double leafSizeX = 1; // 体素大小
    double leafSizeY = 1;
    double leafSizeZ = 1;
    auto vtk_points_denoise = denoisePointCloud(polyData_points,cleanTolerance,leafSizeX,leafSizeY,leafSizeZ);
    return std::make_tuple(vtk_points_denoise,filtered_points);
}