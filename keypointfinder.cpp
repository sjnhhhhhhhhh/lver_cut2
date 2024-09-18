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


#include <iostream>
#include <vector>
#include <cmath>
#include <tuple>

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

#include <pcl/point_types.h>
#include <pcl/filters/voxel_grid.h>

using namespace std;

typedef CGAL::Simple_cartesian<double> K;
typedef K::Point_3 Point_3;
typedef CGAL::Delaunay_triangulation_3<K> Delaunay;
typedef K::Segment_3 Segment_3;
typedef K::Plane_3 Plane_3;

struct PointWithLabel {
    Eigen::Vector3d point;
    int label; // 标识点所属的组：1表示第一组，2表示第二组
};

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
template <typename... PolyDataTypes>
vtkSmartPointer<vtkPolyData> mergePolyData(const PolyDataTypes&... inputs) {
    // 创建一个新的 vtkPoints 用于存储合并后的点
    vtkSmartPointer<vtkPoints> mergedPoints = vtkSmartPointer<vtkPoints>::New();

    // Lambda 表达式用于提取每个 vtkPolyData 的点并添加到 mergedPoints 中
    auto extractAndInsertPoints = [&](const vtkSmartPointer<vtkPolyData>& input) {
        for (vtkIdType i = 0; i < input->GetNumberOfPoints(); ++i) {
            double point[3];
            input->GetPoint(i, point);
            mergedPoints->InsertNextPoint(point);
        }
    };

    // 使用折叠表达式，遍历所有传入的 vtkPolyData 对象
    (extractAndInsertPoints(inputs), ...);

    // 创建一个新的 vtkPolyData 对象并设置合并后的点
    vtkSmartPointer<vtkPolyData> mergedPolyData = vtkSmartPointer<vtkPolyData>::New();
    mergedPolyData->SetPoints(mergedPoints);

    return mergedPolyData;
}

std::vector<Eigen::Vector3d> uniformSampling(const std::vector<Eigen::Vector3d>& points, size_t num_samples) {
    std::vector<Eigen::Vector3d> sampled_points;
    size_t step = points.size() / num_samples;

    for (size_t i = 0; i < points.size(); i += step) {
        sampled_points.push_back(points[i]);
    }

    return sampled_points;
}

// 定义高斯 RBF 核函数
double gaussianRBF(const Eigen::Vector3d& p1, const Eigen::Vector3d& p2, double sigma) {
    return exp(-((p1 - p2).squaredNorm()) / (2 * sigma * sigma));
}

// RBF 拟合函数
Eigen::VectorXd fitRBF(const std::vector<Eigen::Vector3d>& points, double sigma) {
    int n = points.size();
    Eigen::MatrixXd A(n, n);
    Eigen::VectorXd Z(n);

    // 构建矩阵 A 和 z 值
    for (int i = 0; i < n; ++i) {
        Z(i) = points[i][2];  // Z 是 z 值
        for (int j = 0; j < n; ++j) {
            A(i, j) = gaussianRBF(points[i], points[j], sigma);
        }
    }

    // 计算 RBF 权重
    Eigen::VectorXd weights = A.colPivHouseholderQr().solve(Z);

    return weights;  // 返回拟合的 RBF 权重
}

// 使用 RBF 进行插值
double interpolateRBF(const Eigen::VectorXd& weights, const std::vector<Eigen::Vector3d>& points, const Eigen::Vector3d& queryPoint, double sigma) {
    double z = 0.0;
    for (size_t i = 0; i < points.size(); ++i) {
        z += weights(i) * gaussianRBF(queryPoint, points[i], sigma);
    }
    return z;
}

vtkSmartPointer<vtkPolyData> visualizeRBF(const Eigen::VectorXd& weights, const std::vector<Eigen::Vector3d>& points, double bounds[6], int resolution, double sigma) {
    // 1. 创建 VTK 点集
    vtkSmartPointer<vtkPoints> points_grid = vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkCellArray> triangles = vtkSmartPointer<vtkCellArray>::New();

    double xMin = bounds[0];
    double xMax = bounds[1];
    double yMin = bounds[2];
    double yMax = bounds[3];

    // 网格划分
    double dx = (xMax - xMin) / (resolution - 1);
    double dy = (yMax - yMin) / (resolution - 1);

    std::vector<std::vector<vtkIdType>> pointIds(resolution, std::vector<vtkIdType>(resolution));

    // 2. 生成网格点并使用 RBF 插值计算 z 值
    for (int i = 0; i < resolution; ++i) {
        for (int j = 0; j < resolution; ++j) {
            double x = xMin + i * dx;
            double y = yMin + j * dy;
            Eigen::Vector3d queryPoint(x, y, 0.0);  // 查询点，z 值初始为 0
            double z = interpolateRBF(weights, points, queryPoint, sigma);  // 插值得到 z 值
            queryPoint[2] = z;  // 更新 z 值

            // 插入点到 VTK 点集中
            pointIds[i][j] = points_grid->InsertNextPoint(queryPoint[0], queryPoint[1], queryPoint[2]);
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
    polyData->SetPoints(points_grid);
    polyData->SetPolys(triangles);

    return polyData;
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




int main() {
    // 读取两个 STL 文件中的点
    vtkSmartPointer<vtkSTLReader> reader1 = vtkSmartPointer<vtkSTLReader>::New();
    vtkSmartPointer<vtkSTLReader> reader2 = vtkSmartPointer<vtkSTLReader>::New();
    vtkSmartPointer<vtkSTLReader> reader3 = vtkSmartPointer<vtkSTLReader>::New();
    vtkSmartPointer<vtkSTLReader> reader4 = vtkSmartPointer<vtkSTLReader>::New();
    vtkSmartPointer<vtkSTLReader> reader5 = vtkSmartPointer<vtkSTLReader>::New();

    reader1->SetFileName("C:/code/liver_cut/data/rb1v.stl");
    reader2->SetFileName("C:/code/liver_cut/data/rb2v.stl");
    reader3->SetFileName("C:/code/liver_cut/data/rb3v.stl");
    reader4->SetFileName("C:/code/liver_cut/data/rf_up_v.stl");
    reader5->SetFileName("C:/code/liver_cut/data/rf_down_v.stl");
    reader1->Update();
    reader2->Update();
    reader3->Update();
    reader4->Update();
    reader5->Update();

    auto rb1v = reader1->GetOutput();
    auto rb2v = reader2->GetOutput();
    auto rb3v = reader3->GetOutput();
    auto rf_up_v = reader4->GetOutput();
    auto rf_down_v = reader5->GetOutput();
    auto exc_rfd = mergePolyData(rb1v,rb3v,rf_up_v,rf_down_v);

    
    double boundsv1[6];
    double boundsv2[6];
    exc_rfd->GetBounds(boundsv1);
    cout<<"r1v1 range:\n"<<"x:("<<boundsv1[0]<<","<<boundsv1[1]<<")\n";
    cout<<"y:("<<boundsv1[2]<<","<<boundsv1[3]<<")\n";
    cout<<"z:("<<boundsv1[4]<<","<<boundsv1[5]<<")\n";
    rb2v->GetBounds(boundsv2);
    cout<<"r1v1 range:\n"<<"x:("<<boundsv2[0]<<","<<boundsv2[1]<<")\n";
    cout<<"y:("<<boundsv2[2]<<","<<boundsv2[3]<<")\n";
    cout<<"z:("<<boundsv2[4]<<","<<boundsv2[5]<<")\n";

    double mbounds[6];
    max_bounds(boundsv1,boundsv2,mbounds);

    auto pointsvector1 = extractpoints(exc_rfd);
    auto pointsvector2 = extractpoints(rb2v);

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
    points_for_fitting = uniformSampling(points_for_fitting,1000);
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

    double bounds[6];
    polyData_points->GetBounds(bounds);
    cout<<"plane range:\n"<<"x:("<<bounds[0]<<","<<bounds[1]<<")\n";
    cout<<"y:("<<bounds[2]<<","<<bounds[3]<<")\n";

    vtkSmartPointer<vtkPolyDataWriter> writer = vtkSmartPointer<vtkPolyDataWriter>::New();
    writer->SetFileName("C:/code/liver_cut/src/points.vtk");
    writer->SetInputData(polyData_points);
    writer->Write();

    std::cout << "STL file has been written to fitted_surface.stl" << std::endl;

    // 3. 为点创建顶点表示
    vtkSmartPointer<vtkVertexGlyphFilter> glyphFilter = vtkSmartPointer<vtkVertexGlyphFilter>::New();
    glyphFilter->SetInputData(polyData_points);
    glyphFilter->Update();

    // 4. 创建映射器和演员
    vtkSmartPointer<vtkPolyDataMapper> mapper_points = vtkSmartPointer<vtkPolyDataMapper>::New();
    mapper_points->SetInputConnection(glyphFilter->GetOutputPort());

    vtkSmartPointer<vtkActor> actor_points = vtkSmartPointer<vtkActor>::New();
    actor_points->SetMapper(mapper_points);
    actor_points->GetProperty()->SetColor(1.0, 0.0, 0.0);  // 将点集设置为红色
    actor_points->GetProperty()->SetPointSize(2);  // 设置点的大小
/*
    // 体素大小
    double voxelSize = 1.0;
    // 定义体数据网格大小（根据你的点云数据调整大小）
    int gridSize[3] = {100, 100, 100};

    // 1. 将你的 filtered_points 转换为 vtkImageData
    vtkSmartPointer<vtkImageData> imageData = convertPointsToImageData(filtered_points, voxelSize, gridSize);

    // 2. 使用 Marching Cubes 算法提取等值面
    vtkSmartPointer<vtkMarchingCubes> marchingCubes = vtkSmartPointer<vtkMarchingCubes>::New();
    marchingCubes->SetInputData(imageData);
    marchingCubes->SetValue(0, 0.5); // 等值面值为 0.5，通常介于 0 和 1 之间
    marchingCubes->Update();

    // 获取生成的等值面
    vtkSmartPointer<vtkPolyData> polyData = marchingCubes->GetOutput();

    // 创建一个 PolyDataMapper 和 Actor
    vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    mapper->SetInputData(polyData);

    vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
    actor->SetMapper(mapper);
    actor->GetProperty()->SetOpacity(0.7);
*/
    
    vtkSmartPointer<vtkSurfaceReconstructionFilter> surfaceReconstruction = vtkSmartPointer<vtkSurfaceReconstructionFilter>::New();
    surfaceReconstruction->SetInputData(polyData_points);
    surfaceReconstruction->Update();

    vtkSmartPointer<vtkContourFilter> contourFilter = vtkSmartPointer<vtkContourFilter>::New();
    contourFilter->SetInputConnection(surfaceReconstruction->GetOutputPort());
    contourFilter->SetValue(0, 0.0);
    contourFilter->Update();

    vtkSmartPointer<vtkPolyDataMapper> mapper_vtk = vtkSmartPointer<vtkPolyDataMapper>::New();
    mapper_vtk->SetInputConnection(contourFilter->GetOutputPort());
    vtkSmartPointer<vtkActor> actor_vtk = vtkSmartPointer<vtkActor>::New();
    actor_vtk->SetMapper(mapper_vtk);


    // 创建第一个映射器和演员
    vtkSmartPointer<vtkPolyDataMapper> mapper1 = vtkSmartPointer<vtkPolyDataMapper>::New();
    mapper1->SetInputConnection(reader1->GetOutputPort());

    vtkSmartPointer<vtkActor> actor1 = vtkSmartPointer<vtkActor>::New();
    actor1->SetMapper(mapper1);
    actor1->GetProperty()->SetOpacity(0.7);  // 设置透明度

    // 创建第二个映射器和演员
    vtkSmartPointer<vtkPolyDataMapper> mapper2 = vtkSmartPointer<vtkPolyDataMapper>::New();
    mapper2->SetInputConnection(reader2->GetOutputPort());

    vtkSmartPointer<vtkActor> actor2 = vtkSmartPointer<vtkActor>::New();
    actor2->SetMapper(mapper2);
    actor2->GetProperty()->SetOpacity(0.7);  // 设置透明度

    // 创建第二个映射器和演员
    vtkSmartPointer<vtkPolyDataMapper> mapper3 = vtkSmartPointer<vtkPolyDataMapper>::New();
    mapper3->SetInputConnection(reader3->GetOutputPort());

    vtkSmartPointer<vtkActor> actor3 = vtkSmartPointer<vtkActor>::New();
    actor3->SetMapper(mapper3);
    actor3->GetProperty()->SetOpacity(0.7);  // 设置透明度

    // 创建第二个映射器和演员
    vtkSmartPointer<vtkPolyDataMapper> mapper4 = vtkSmartPointer<vtkPolyDataMapper>::New();
    mapper4->SetInputConnection(reader4->GetOutputPort());

    vtkSmartPointer<vtkActor> actor4 = vtkSmartPointer<vtkActor>::New();
    actor4->SetMapper(mapper4);
    actor4->GetProperty()->SetOpacity(0.7);  // 设置透明度






    // 创建第二个映射器和演员
    vtkSmartPointer<vtkPolyDataMapper> mapper5 = vtkSmartPointer<vtkPolyDataMapper>::New();
    mapper5->SetInputConnection(reader5->GetOutputPort());

    vtkSmartPointer<vtkActor> actor5 = vtkSmartPointer<vtkActor>::New();
    actor5->SetMapper(mapper5);
    actor5->GetProperty()->SetOpacity(0.7);  // 设置透明度

    // 创建一个平面源，并设置其方向和位置
    vtkSmartPointer<vtkPlaneSource> planeSource = vtkSmartPointer<vtkPlaneSource>::New();
    planeSource->SetCenter(centroid[0], centroid[1], centroid[2]);
    planeSource->SetNormal(normal[0], normal[1], normal[2]);

    
    planeSource->Update();


    // 创建渲染器、渲染窗口和交互器
    vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
    vtkSmartPointer<vtkRenderWindow> renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
    renderWindow->SetWindowName("STL Model and Plane Display");
    renderWindow->AddRenderer(renderer);

    vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
    renderWindowInteractor->SetRenderWindow(renderWindow);

    // 将演员添加到场景中
    renderer->AddActor(actor1);  // 第一个 STL 模型
    renderer->AddActor(actor2);  // 第二个 STL 模型
    renderer->AddActor(actor3);
    renderer->AddActor(actor4);
    renderer->AddActor(actor5);
    //renderer->AddActor(actor);
    renderer->AddActor(actor_vtk);  // 拟合平面
    renderer->AddActor(actor_points);
    //renderer->AddActor(actor_rbf);

    //renderer->AddActor(rbf_actor);

    // 设置背景颜色
    renderer->SetBackground(0.1, 0.2, 0.3);

    // 调整摄像机视角，使得模型和坐标轴都能正确显示
    renderer->ResetCamera();
    renderWindow->Render();
    renderWindowInteractor->Start();

    return EXIT_SUCCESS;
}

    /*
    目前尝试了四次拟合，仍然不理想，接下来可以尝试下vtk的行进立方，或者其他什么办法

    或者换一下切割顺序，简化一下右前下段的切割面复杂度
    */    

    /*auto polydata = vec_to_poly(points_for_fitting);
    double bounds[6];
    polydata->GetBounds(bounds);
    cout<<"plane range:\n"<<"x:("<<bounds[0]<<","<<bounds[1]<<")\n";
    cout<<"y:("<<bounds[2]<<","<<bounds[3]<<")\n";

——————————————————————————稠密点云泊松重建复杂曲面——————————————————————————————————————
    /*void reconstructSurface(const pcl::PointCloud<pcl::PointXYZ>::Ptr& cloud) {
    // 1. 创建 Poisson 重建对象
    pcl::Poisson<pcl::PointXYZ> poisson;
    poisson.setInputCloud(cloud);

    // 2. 执行重建，生成 PolygonMesh
    pcl::PolygonMesh mesh;
    poisson.reconstruct(mesh);

    // 3. 保存生成的表面为 PLY 文件
    pcl::io::savePLYFile("C:/code/liver_cut/src/reconstructed_surface.ply", mesh);

    std::cout << "Reconstructed surface saved to reconstructed_surface.ply" << std::endl;
}*/

// 将 std::vector<Eigen::Vector3d> 转换为 pcl::PointCloud<pcl::PointXYZ>::Ptr
/*pcl::PointCloud<pcl::PointXYZ>::Ptr convertToPCLPointCloud(const std::vector<Eigen::Vector3d>& eigen_points) {
    pcl::PointCloud<pcl::PointXYZ>::Ptr cloud(new pcl::PointCloud<pcl::PointXYZ>);

    for (const auto& point : eigen_points) {
        pcl::PointXYZ pcl_point(point[0], point[1], point[2]);
        cloud->points.push_back(pcl_point);
    }

    cloud->width = cloud->points.size();
    cloud->height = 1;
    cloud->is_dense = true;

    return cloud;
}*/
    //auto cloud = convertToPCLPointCloud(filtered_points);
    //reconstructSurface(cloud);

    //——————————————————————拟合平面——————————————————————————————


    /*auto beta = fitPolynomialSurface(filtered_points);
    auto polyData_plane = visualizePolynomialSurface(beta,mbounds,50);

    vtkSmartPointer<vtkPolyDataMapper> mapper_plane = vtkSmartPointer<vtkPolyDataMapper>::New();
    mapper_plane->SetInputData(polyData_plane);

    vtkSmartPointer<vtkActor> actor_plane = vtkSmartPointer<vtkActor>::New();
    actor_plane->SetMapper(mapper_plane);
    actor_plane->GetProperty()->SetOpacity(0.7);  // 设置透明度
    actor_plane->GetProperty()->SetColor(1.0,0.0,0.0);  // 设置透明度*/



/*
    // 使用 Delaunay3D 进行三角剖分
    vtkSmartPointer<vtkDelaunay3D> delaunay = vtkSmartPointer<vtkDelaunay3D>::New();
    delaunay->SetInputData(polyData);
    delaunay->Update();

    // 4. 使用 vtkGeometryFilter 将 UnstructuredGrid 转换为 PolyData
    vtkSmartPointer<vtkGeometryFilter> geometryFilter = vtkSmartPointer<vtkGeometryFilter>::New();
    geometryFilter->SetInputData(delaunay->GetOutput());  // 输入 Delaunay 的输出 (UnstructuredGrid)
    geometryFilter->Update();  // 执行转换
    auto beta = fitPolynomialSurface(points_for_fitting);
    auto polydata2 = visualizePolynomialSurface(beta,bounds,50);
    //5. 使用 vtkSTLWriter 将 PolyData 写入 STL 文件
    vtkSmartPointer<vtkSTLWriter> writer = vtkSmartPointer<vtkSTLWriter>::New();
    writer->SetFileName("C:/code/liver_cut/src/fitted_surface.stl");  // 输出的 STL 文件路径和文件名
    writer->SetInputData(polydata2);  // 将转换后的 PolyData 作为输入
    writer->Write();  // 执行写入

    std::cout << "STL file has been written to fitted_surface.stl" << std::endl;

    // 创建映射器和演员用于可视化
    vtkSmartPointer<vtkDataSetMapper> mapper_del = vtkSmartPointer<vtkDataSetMapper>::New();
    mapper_del->SetInputData(polydata2);

    vtkSmartPointer<vtkActor> actor_del = vtkSmartPointer<vtkActor>::New();
    actor_del->SetMapper(mapper_del);
    actor_del->GetProperty()->SetColor(1.0, 0.5, 0.0);  // 设置表面颜色
    actor_del->GetProperty()->SetOpacity(0.7);  // 设置透明度*/

    /* 创建平面映射器和演员
    vtkSmartPointer<vtkPolyDataMapper> planeMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    planeMapper->SetInputConnection(planeSource->GetOutputPort());

    vtkSmartPointer<vtkActor> planeActor = vtkSmartPointer<vtkActor>::New();
    planeActor->SetMapper(planeMapper);
    planeActor->GetProperty()->SetColor(1.0, 0.0, 0.0);  // 设置平面为红色
    planeActor->GetProperty()->SetOpacity(0.5);  // 设置透明度*/

    // 平面拟合
    /*if (points_for_fitting.size() > 3) {
        // 计算点集的质心
        
        for (const auto& point : points_for_fitting) {
            centroid += point;
        }
        centroid /= points_for_fitting.size();

        // 计算协方差矩阵
        Eigen::Matrix3d covariance_matrix = Eigen::Matrix3d::Zero();
        for (const auto& point : points_for_fitting) {
            Eigen::Vector3d centered_point = point - centroid;
            covariance_matrix += centered_point * centered_point.transpose();
        }

        // 使用特征值分解计算协方差矩阵的特征向量
        Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eigen_solver(covariance_matrix);
        normal = eigen_solver.eigenvectors().col(0);  // 最小特征值对应的特征向量是平面的法向量

        std::cout << "Fitted plane normal: " << normal.transpose() << std::endl;
        std::cout << "Plane passes through point: " << centroid.transpose() << std::endl;
    } else {
        std::cout << "Not enough points to fit a plane." << std::endl;
    }*/
/*整体思路：
    先提取stl点，然后全部合并三角化，然后找到两端标签不同的线段，这些就是边界了，然后提取这些线段的voronoi对偶线段，它们端点组成的点云就是我们要找的分割平面的点集
    
    -->随后我们只需要进行稀疏化，然后拟合曲面就可以了
    
    在拟合玩所以曲面后，我们可以按照顺序，逐个应用曲面给cutter，然后切出来成果

*/

//方法一：先找到距离最近的两点，做连线，得到法向量
    //接下来让所有点垂直这个法向量去寻找与另一个点集的交点，算出中点
    //由所有中点拟合一个分割平面
    /*
    现有问题：
    （1）过穿问题：由于我们只想要最靠近的点（即靠左的血管的右面和靠右的血管的左面）之间的中点，怎么避免遍历的时候出现左侧点与另一点集左侧点的计算
    （2）交点问题：还是存在这个问题，如果你用到了交点计算并且是双层嵌套循环，那么时间复杂度是O(n^2)，太蠢了太慢了
    */

//方法二：voronoi图
/*
由于voronoi图天生就是对点的区域划分，所以如果我们能生成维诺图那么平面也可以找到了
因为我们可以生成连续的分割线段，然后就可以得到连续线段
*//*int main(){
    // 提取两个STL文件的点集
    vtkSmartPointer<vtkSTLReader> reader1 = vtkSmartPointer<vtkSTLReader>::New();
    vtkSmartPointer<vtkSTLReader> reader2 = vtkSmartPointer<vtkSTLReader>::New();
    reader1->SetFileName("C:/code/liver_cut/data/rb1v.stl");
    reader2->SetFileName("C:/code/liver_cut/data/rb2v.stl");
    reader1->Update();
    reader2->Update();

    auto input1 = reader1->GetOutput();
    auto input2 = reader2->GetOutput();

    auto pointsvector1 = extractpoints(input1);
    auto pointsvector2 = extractpoints(input2);

    cout << "number of rb1v points: " << pointsvector1.size() << endl;
    cout << "number of rb2v points: " << pointsvector2.size() << "\n";

    // 合并点集并标记点的归属
    auto points = mergePointSets(pointsvector1, pointsvector2);

    // 创建 CGAL 的 Point_3 类型点集
    std::vector<Point_3> cgal_points;
    for (const auto& pointWithLabel : points) {
        cgal_points.push_back(Point_3(pointWithLabel.point[0], pointWithLabel.point[1], pointWithLabel.point[2]));
    }

    cout<<"generating Delaunay------\n";


    // 生成 Delaunay 三角化
    Delaunay dt;
    dt.insert(cgal_points.begin(), cgal_points.end());

    // 提取分割线段
    std::vector<std::pair<Point_3, Point_3>> splitting_lines;

    cout<<"extracting splitting lines\n";
    // 遍历 Delaunay 三角化中的每条边，判断其是否为分割边界
    for (auto it = dt.finite_edges_begin(); it != dt.finite_edges_end(); ++it) {
        auto cell = it->first;
        int i = it->second;  // 第一个顶点的索引
        int j = it->third;   // 第二个顶点的索引

        auto vertex1 = cell->vertex(i);
        auto vertex2 = cell->vertex(j);

        // 查找点在原始点集中的位置并获取标签
        int index1 = std::distance(cgal_points.begin(), std::find(cgal_points.begin(), cgal_points.end(), vertex1->point()));
        int index2 = std::distance(cgal_points.begin(), std::find(cgal_points.begin(), cgal_points.end(), vertex2->point()));
        int label1 = points[index1].label;
        int label2 = points[index2].label;

        // 如果边的两个顶点分别来自不同的点集
        if (label1 != label2) {
            Point_3 p1 = vertex1->point();
            Point_3 p2 = vertex2->point();

            // 将分割线段加入分割线段集合
            splitting_lines.push_back({p1, p2});
        }
    }

    cout << "Number of splitting lines: " << splitting_lines.size() << endl;

    // Voronoi 图提取：提取 Voronoi 图中与 Delaunay 边对偶的结构
    for (auto it = dt.finite_facets_begin(); it != dt.finite_facets_end(); ++it) {
        CGAL::Object voronoi_object = dt.dual(*it); // Voronoi 对偶（可能是平面、点等）

        //dual函数有多种重载，传入面返回对偶面，传入线段返回对偶线段



        // 尝试将对偶对象转换为线段
        Segment_3 voronoi_segment;
        if (CGAL::assign(voronoi_segment, voronoi_object)) {
            cout << "Voronoi segment (dual): " << voronoi_segment << endl;
        } else {
            cout << "Voronoi dual object is not a segment." << endl;
        }
    }

    return 0;
    }*/