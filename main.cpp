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
#include <vtkSurfaceReconstructionFilter.h>
#include <vtkContourFilter.h>
#include <vtkAppendPolyData.h>
#include <vtkImplicitPolyDataDistance.h>
#include <vtkClipPolyData.h>

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

// 拟合平面
vtkSmartPointer<vtkImplicitPolyDataDistance> plane_generator(std::vector<Eigen::Vector3d> filtered_points){
    vtkSmartPointer<vtkPoints> vtk_points = vtkSmartPointer<vtkPoints>::New();
    for (const auto& point : filtered_points) {
        vtk_points->InsertNextPoint(point[0], point[1], point[2]);
    }

    vtkSmartPointer<vtkPolyData> polyData_points = vtkSmartPointer<vtkPolyData>::New();
    polyData_points->SetPoints(vtk_points);

    vtkSmartPointer<vtkSurfaceReconstructionFilter> surfaceReconstruction = vtkSmartPointer<vtkSurfaceReconstructionFilter>::New();
    surfaceReconstruction->SetInputData(polyData_points);
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




// 切割顺序：尾状叶段-->左外叶上段-->左外叶下段-->左内叶段-->

int main(){
    vtkSmartPointer<vtkSTLReader> reader1 = vtkSmartPointer<vtkSTLReader>::New();
    reader1->SetFileName("C:/code/liver_cut/data/rb1v.stl");
    reader1->Update();
    auto rb1v = reader1->GetOutput();
    std::cout<<"rb1v:"<<rb1v->GetNumberOfPoints()<<"\n";

    vtkSmartPointer<vtkSTLReader> reader2 = vtkSmartPointer<vtkSTLReader>::New();
    reader2->SetFileName("C:/code/liver_cut/data/rb2v.stl");
    reader2->Update();
    auto rb2v = reader2->GetOutput();
    std::cout<<rb2v->GetNumberOfPoints()<<"\n";

    vtkSmartPointer<vtkSTLReader> reader3 = vtkSmartPointer<vtkSTLReader>::New();
    reader3->SetFileName("C:/code/liver_cut/data/rb3v.stl");
    reader3->Update();
    auto rb3v = reader3->GetOutput();
    
    vtkSmartPointer<vtkSTLReader> reader4 = vtkSmartPointer<vtkSTLReader>::New();
    reader4->SetFileName("C:/code/liver_cut/data/rf_up_v.stl");
    reader4->Update();
    auto rf_up_v = reader4->GetOutput();

    vtkSmartPointer<vtkSTLReader> reader5 = vtkSmartPointer<vtkSTLReader>::New();
    reader5->SetFileName("C:/code/liver_cut/data/rf_down_v.stl");
    reader5->Update();
    auto rf_down_v = reader5->GetOutput();

    vtkSmartPointer<vtkSTLReader> reader6 = vtkSmartPointer<vtkSTLReader>::New();
    reader6->SetFileName("C:/code/liver_cut/data/left_inside_v.stl");
    reader6->Update();
    auto left_inside_v = reader6->GetOutput();

    vtkSmartPointer<vtkSTLReader> reader7 = vtkSmartPointer<vtkSTLReader>::New();
    reader7->SetFileName("C:/code/liver_cut/data/left_ou_v.stl");
    reader7->Update();
    auto left_ou_v = reader7->GetOutput();

    vtkSmartPointer<vtkSTLReader> reader8 = vtkSmartPointer<vtkSTLReader>::New();
    reader8->SetFileName("C:/code/liver_cut/data/left_od_v.stl");
    reader8->Update();
    auto left_od_v = reader8->GetOutput();

    vtkSmartPointer<vtkSTLReader> reader9 = vtkSmartPointer<vtkSTLReader>::New();
    reader9->SetFileName("C:/code/liver_cut/data/tail_v.stl");
    reader9->Update();
    auto tail = reader9->GetOutput();

    // 由于没给我整个肝脏的stl，我只能自己合并下了，你要是有的话就删了吧
    vtkSmartPointer<vtkPolyData> liver = vtkSmartPointer<vtkPolyData>::New();
        vtkSmartPointer<vtkSTLReader> reader10 = vtkSmartPointer<vtkSTLReader>::New();
        reader10->SetFileName("C:/code/liver_cut/data/1.stl");
        reader10->Update();
        liver = mergePolyData(liver,reader10->GetOutput());

        reader10->SetFileName("C:/code/liver_cut/data/2.stl");
        reader10->Update();
        liver = mergePolyData(liver,reader10->GetOutput());

        reader10->SetFileName("C:/code/liver_cut/data/3.stl");
        reader10->Update();
        liver = mergePolyData(liver,reader10->GetOutput());

        reader10->SetFileName("C:/code/liver_cut/data/4a.stl");
        reader10->Update();
        liver = mergePolyData(liver,reader10->GetOutput());

        reader10->SetFileName("C:/code/liver_cut/data/4b.stl");
        reader10->Update();
        liver = mergePolyData(liver,reader10->GetOutput());

        reader10->SetFileName("C:/code/liver_cut/data/5.stl");
        reader10->Update();
        liver = mergePolyData(liver,reader10->GetOutput());

        reader10->SetFileName("C:/code/liver_cut/data/6.stl");
        reader10->Update();
        liver = mergePolyData(liver,reader10->GetOutput());
        
        reader10->SetFileName("C:/code/liver_cut/data/7.stl");
        reader10->Update();
        liver = mergePolyData(liver,reader10->GetOutput());

        reader10->SetFileName("C:/code/liver_cut/data/8.stl");
        reader10->Update();
        liver = mergePolyData(liver,reader10->GetOutput());   

        std::cout<<liver->GetNumberOfPoints()<<"\n";


    // 第一刀，尾状段
    auto exc_tail = mergePolyData(rb1v,rb2v,rb3v,rf_up_v,rf_down_v,left_inside_v,left_ou_v,left_od_v);
    auto [cut_tail,cut_tail_points] = find_points(exc_tail,tail);
    auto imp1 = plane_generator(cut_tail_points);
    auto [liver1,tail_poly] = cut_run2(imp1,liver);
    vtkSmartPointer<vtkSTLWriter> stlWriter = vtkSmartPointer<vtkSTLWriter>::New();
    stlWriter->SetFileName("C:/code/liver_cut/new_stl/tail.stl");
    stlWriter->SetInputData(tail_poly);
    stlWriter->Write();
    std::cout << "Successfully saved" << std::endl;
    std::cout<<"\ntail cut completed\n";

    // 第二刀，左外叶上段
    auto exc_left_ou_v = mergePolyData(left_inside_v,left_od_v,tail);
    auto [cut_left_ou,cut_left_ou_points] = find_points(exc_left_ou_v,left_ou_v); 
    auto imp2 = plane_generator(cut_left_ou_points);
    auto [liver2,left_ou_poly] = cut_run2(imp2,liver1);
    vtkSmartPointer<vtkSTLWriter> stlWriter2 = vtkSmartPointer<vtkSTLWriter>::New();
    stlWriter2->SetFileName("C:/code/liver_cut/new_stl/left_ou.stl");
    stlWriter2->SetInputData(left_ou_poly);
    stlWriter2->Write();
    std::cout << "Successfully saved" << std::endl;
    std::cout<<"left_ou cut completed\n";

    // 第三刀，左外叶下段
    auto [cut_left_od,cut_left_od_points] = find_points(left_od_v,left_inside_v); 
    auto imp3 = plane_generator(cut_left_od_points);
    auto [liver3,left_od_poly] = cut_run2(imp3,liver2);

    std::cout<<"\nleftod cut completed\n";

    // 第四刀，左内叶段
    auto exc_left_inside_v = mergePolyData(rf_up_v,rf_down_v);
    auto [cut_left_inside,cut_left_inside_points] = find_points(exc_left_inside_v,left_inside_v); 
    auto imp4 = plane_generator(cut_left_inside_points);
    auto [liver4,left_inside_poly] = cut_run2(imp4,liver3);

    std::cout<<"\nleft_inside cut completed\n";

    // 第五刀，右后3段
    auto exc_rb3v = mergePolyData(rb1v,rb2v,rf_down_v,rf_up_v);
    auto [rb3b,rb3b_points] = find_points(exc_rb3v,rb3v);
    auto imp5 = plane_generator(rb3b_points);
    auto [liver5,rb3b_poly] = cut_run2(imp5,liver4);

    std::cout<<"\nrb3b cut completed\n";

    // 第六刀，右后2段
    auto exc_rb2v = mergePolyData(rb1v,rf_up_v);
    auto [rb2b,rb2b_points] = find_points(exc_rb2v,rb2v);
    auto imp6 = plane_generator(rb2b_points);
    auto [liver6,rb2b_poly] = cut_run2(imp6,liver5);

    std::cout<<"\nrb2b cut completed\n";

    // 第七刀，右后1段
    auto exc_rb1v = mergePolyData(rf_down_v,rf_up_v);
    auto [rb1b,rb1b_points] = find_points(exc_rb1v,rb1v);
    auto imp7 = plane_generator(rb3b_points);
    auto [liver7,rb1b_poly] = cut_run2(imp7,liver6);

    std::cout<<"\nrb1b cut completed\n";

    // 第八刀，右前上段
    auto exc_rf_up_v = rf_down_v;
    auto [rf_up_b,rf_up_v_points] = find_points(exc_rf_up_v,rf_up_v);
    auto imp8 = plane_generator(rf_up_v_points);
    auto [liver8,rf_up_v_poly] = cut_run2(imp8,liver7);

    std::cout<<"\nrf_up cut completed\n";

    // 最后剩下右前下段
    auto rf_down_v_poly = liver8;

    std::cout<<"\nrf_down cut completed\n";

    return 0;




}
/*
    // 第一刀，尾状段
    auto exc_tail = mergePolyData(rb1v,rb2v,rb3v,rf_up_v,rf_down_v,left_inside_v,left_ou_v,left_od_v);
    auto [cut_tail,cut_tail_points] = find_points(exc_tail,tail);
    auto imp1 = plane_generator(cut_tail_points);
    auto tail_poly = cut_run(imp1,liver);

    // 第二刀，左后上段
    auto exc_left_ou_v = mergePolyData(left_inside_v,left_od_v,tail);
    auto [cut_left_ou,cut_left_ou_points] = find_points(exc_left_ou_v,left_ou_v); 
    auto imp2 = plane_generator(cut_left_ou_points);
    auto left_ou_poly = cut_run(imp2,liver);

    // 第三刀，左后下段
    auto [cut_left_od,cut_left_od_points] = find_points(left_od_v,left_inside_v); 
    auto imp3 = plane_generator(cut_left_od_points);
    auto left_od_poly = cut_run(imp3,liver);

    // 第四刀，左内叶段
    auto exc_left_inside_v = mergePolyData(rf_up_v,rf_down_v);
    auto [cut_left_inside,cut_left_inside_points] = find_points(exc_left_inside_v,left_inside_v); 
    auto imp4 = plane_generator(cut_left_inside_points);
    auto left_inside_poly = cut_run(imp4,liver);

    // 第五刀，右后3段
    auto exc_rb3v = mergePolyData(rb1v,rb2v,rf_down_v,rf_up_v,tail,left_inside_v);
    auto [rb3b,rb3b_points] = find_points(exc_rb3v,rb3v);
    auto imp5 = plane_generator(rb3b_points);
    auto rb3b_poly = cut_run(imp5,liver);

    // 第六刀，右后3段
    auto exc_rb3v = mergePolyData(rb1v,rb2v,rf_down_v,rf_up_v,tail,left_inside_v);
    auto [rb3b,rb3b_points] = find_points(exc_rb3v,rb3v);
    auto imp5 = plane_generator(rb3b_points);
    auto rb3b_poly = cut_run(imp5,liver);

    // 第五刀，右后3段
    auto exc_rb3v = mergePolyData(rb1v,rb2v,rf_down_v,rf_up_v,tail,left_inside_v);
    auto [rb3b,rb3b_points] = find_points(exc_rb3v,rb3v);
    auto imp5 = plane_generator(rb3b_points);
    auto rb3b_poly = cut_run(imp5,liver);

    // 第五刀，右后3段
    auto exc_rb3v = mergePolyData(rb1v,rb2v,rf_down_v,rf_up_v,tail,left_inside_v);
    auto [rb3b,rb3b_points] = find_points(exc_rb3v,rb3v);
    auto imp5 = plane_generator(rb3b_points);
    auto rb3b_poly = cut_run(imp5,liver);
*/
