//根据stl集，写一下自动化的合并函数，然后依次应用find_points

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
#include <utility>  // For std::move

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


// 可变参数模板函数，直接完成 vtkPolyData 合并并转换为 Eigen::Vector3d
template <typename... PolyDataTypes>
std::vector<Eigen::Vector3d> mergeAndExtractPoints(const PolyDataTypes&... inputs) {
    std::vector<Eigen::Vector3d> points;

    // Lambda 表达式提取点集并转换为 Eigen::Vector3d，直接加入 points
    auto extractAndInsertPoints = [&](const vtkSmartPointer<vtkPolyData>& input) {
        for (vtkIdType i = 0; i < input->GetNumberOfPoints(); ++i) {
            double point[3];
            input->GetPoint(i, point);
            points.emplace_back(Eigen::Vector3d(point[0], point[1], point[2]));
        }
    };

    // 使用折叠表达式遍历所有传入的 vtkPolyData 对象
    (extractAndInsertPoints(inputs), ...);

    return points;
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

std::vector<Eigen::Vector3d> extract_points_from_poly(vtkSmartPointer<vtkPolyData> input){
    std::vector<Eigen::Vector3d> points;
    for(vtkIdType i = 0;i < input->GetNumberOfPoints();++i){
        double point[3];
        input->GetPoint(i,point);
        Eigen::Vector3d point_vec(point[0],point[1],point[2]);
        points.push_back(point_vec);
    }
    return points;
}

// 切割顺序：尾状叶段-->左外叶上段-->左外叶下段-->左内叶段-->

int main(){
    vtkSmartPointer<vtkSTLReader> reader1 = vtkSmartPointer<vtkSTLReader>::New();
    reader1->SetFileName("C:/code/liver_cut/data/rb1v.stl");
    auto rb1v = reader1->GetOutput();

    vtkSmartPointer<vtkSTLReader> reader2 = vtkSmartPointer<vtkSTLReader>::New();
    reader2->SetFileName("C:/code/liver_cut/data/rb2v.stl");
    auto rb2v = reader2->GetOutput();

    vtkSmartPointer<vtkSTLReader> reader3 = vtkSmartPointer<vtkSTLReader>::New();
    reader3->SetFileName("C:/code/liver_cut/data/rb3v.stl");
    auto rb3v = reader3->GetOutput();
    
    vtkSmartPointer<vtkSTLReader> reader4 = vtkSmartPointer<vtkSTLReader>::New();
    reader4->SetFileName("C:/code/liver_cut/data/rf_up_v.stl");
    auto rf_up_v = reader4->GetOutput();

    vtkSmartPointer<vtkSTLReader> reader5 = vtkSmartPointer<vtkSTLReader>::New();
    reader5->SetFileName("C:/code/liver_cut/data/rf_down_v.stl");
    auto rf_down_v = reader5->GetOutput();

    vtkSmartPointer<vtkSTLReader> reader6 = vtkSmartPointer<vtkSTLReader>::New();
    reader6->SetFileName("C:/code/liver_cut/data/left_inside_v.stl");
    auto left_inside_v = reader6->GetOutput();

    vtkSmartPointer<vtkSTLReader> reader7 = vtkSmartPointer<vtkSTLReader>::New();
    reader7->SetFileName("C:/code/liver_cut/data/left_ou_v.stl");
    auto left_ou_v = reader7->GetOutput();

    vtkSmartPointer<vtkSTLReader> reader8 = vtkSmartPointer<vtkSTLReader>::New();
    reader8->SetFileName("C:/code/liver_cut/data/left_od_v.stl");
    auto left_od_v = reader8->GetOutput();

    vtkSmartPointer<vtkSTLReader> reader9 = vtkSmartPointer<vtkSTLReader>::New();
    reader9->SetFileName("C:/code/liver_cut/data/tail_v.stl");
    auto tail = reader9->GetOutput();

    // 第一刀，尾状段
    auto exc_tail = mergePolyData(rb1v,rb2v,rb3v,rf_up_v,rf_down_v,left_inside_v,left_ou_v,left_od_v);
    auto [cut_tail,cut_tail_points] = find_points(exc_tail,tail); 

    // 第二刀，左后上段
    auto exc_left_ou_v = mergePolyData(left_inside_v,left_od_v,tail);
    auto [cut_left_ou,cut_left_ou_points] = find_points(exc_left_ou_v,left_ou_v); 

    // 第三刀，左后下段
    auto [cut_left_od,cut_left_od_points] = find_points(left_od_v,left_inside_v); 

    // 第四刀，左内叶段
    auto exc_left_inside_v = mergePolyData(rf_up_v,rf_down_v);
    auto [cut_left_ou,cut_left_ou_points] = find_points(exc_left_inside_v,left_inside_v); 








}