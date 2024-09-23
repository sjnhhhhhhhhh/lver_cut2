// 根据stl集，写一下自动化的合并函数，然后依次应用find_points
// 目前看来切割有问题，平面拟合没问题，只可能是隐含函数和切割出问题了
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
#include <vtkSampleFunction.h>

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


    // 第1刀，尾状段
    auto exc_tail = mergePolyData(rb1v,rb2v,rb3v,rf_up_v,rf_down_v,left_inside_v,left_ou_v,left_od_v);
    auto [cut_tail,cut_tail_points] = find_points(exc_tail,tail);
    auto imp1 = plane_generator(cut_tail);
    auto [liver1,tail_poly] = cut_run2(imp1,liver);
    vtkSmartPointer<vtkSTLWriter> stlWriter = vtkSmartPointer<vtkSTLWriter>::New();
    stlWriter->SetFileName("C:/code/liver_cut/new_stl/tail.stl");
    stlWriter->SetInputData(tail_poly);
    stlWriter->Write();
    std::cout << "Successfully saved" << std::endl;
    std::cout<<"\ntail cut completed\n";

    // 第2刀，右后3段
    auto exc_rb3v = mergePolyData(rb1v,rb2v,rf_down_v,rf_up_v,tail);
    auto [rb3b,rb3b_points] = find_points(exc_rb3v,rb3v);
    auto imp2 = plane_generator(rb3b);
    auto [liver2,rb3b_poly] = cut_run2(imp2,liver1);
    stlWriter->SetFileName("C:/code/liver_cut/new_stl/rb3b.stl");
    stlWriter->SetInputData(rb3b_poly);
    stlWriter->Write();
    std::cout<<"\nrb3b cut completed\n";
    
    // 第3刀，右后2段
    auto exc_rb2v = mergePolyData(rb1v,rb3v,rf_up_v,tail);
    auto [rb2b,rb2b_points] = find_points(exc_rb2v,rb2v);
    auto imp3 = plane_generator(rb2b);
    auto [liver3,rb2b_poly] = cut_run2(imp3,liver2);
    stlWriter->SetFileName("C:/code/liver_cut/new_stl/rb2b.stl");
    stlWriter->SetInputData(rb2b_poly);
    stlWriter->Write();
    std::cout<<"\nrb2b cut completed\n";

    // 第4刀，右后1段
    auto exc_rb1v = mergePolyData(rf_down_v,rf_up_v,rb3v,tail);
    auto [rb1b,rb1b_points] = find_points(exc_rb1v,rb1v);
    auto imp4 = plane_generator(rb3b);
    auto [liver4,rb1b_poly] = cut_run2(imp4,liver3);
    stlWriter->SetFileName("C:/code/liver_cut/new_stl/rb1b.stl");
    stlWriter->SetInputData(rb1b_poly);
    stlWriter->Write();
    std::cout<<"\nrb1b cut completed\n";

    // 第5刀，右前上段
    auto exc_rf_up_v = rf_down_v;
    auto [rf_up_b,rf_up_v_points] = find_points(exc_rf_up_v,rf_up_v);
    auto imp5 = plane_generator(rf_up_v);
    auto [liver5,rf_up_v_poly] = cut_run2(imp5,liver4);
    stlWriter->SetFileName("C:/code/liver_cut/new_stl/rf_up.stl");
    stlWriter->SetInputData(rf_up_v_poly);
    stlWriter->Write();
    std::cout<<"\nrf_up cut completed\n";

    // 第6刀，右前下段
    auto exc_rf_down_v = left_inside_v;
    auto [rf_down_b,rf_down_v_points] = find_points(exc_rf_down_v,rf_down_v);
    auto imp6 = plane_generator(rf_up_v);
    auto [liver6,rf_down_v_poly] = cut_run2(imp6,liver5);
    stlWriter->SetFileName("C:/code/liver_cut/new_stl/rf_down.stl");
    stlWriter->SetInputData(rf_down_v_poly);
    stlWriter->Write();
    std::cout<<"\nrf_down cut completed\n";

    // 第7刀，左内叶段
    auto exc_left_inside_v = mergePolyData(left_ou_v,left_od_v);
    auto [cut_left_inside,cut_left_inside_points] = find_points(exc_left_inside_v,left_inside_v); 
    auto imp7 = plane_generator(cut_left_inside);
    auto [liver7,left_inside_poly] = cut_run2(imp7,liver6);
    stlWriter->SetFileName("C:/code/liver_cut/new_stl/left_inside.stl");
    stlWriter->SetInputData(left_inside_poly);
    stlWriter->Write();
    std::cout << "Successfully saved" << std::endl;
    std::cout<<"\nleft_inside cut completed\n";

    // 第8刀，左外叶上段
    auto exc_left_ou_v = mergePolyData(left_od_v);
    auto [cut_left_ou,cut_left_ou_points] = find_points(exc_left_ou_v,left_ou_v); 
    auto imp8 = plane_generator(cut_left_ou);
    auto [liver8,left_ou_poly] = cut_run2(imp8,liver7);
    vtkSmartPointer<vtkSTLWriter> stlWriter2 = vtkSmartPointer<vtkSTLWriter>::New();
    stlWriter2->SetFileName("C:/code/liver_cut/new_stl/left_ou.stl");
    stlWriter2->SetInputData(left_ou_poly);
    stlWriter2->Write();
    std::cout << "Successfully saved" << std::endl;
    std::cout<<"left_ou cut completed\n";
/*
    // 第三刀，左外叶下段
    auto exc_left_od_v = mergePolyData(left_inside_v,left_ou_v,tail);
    auto [cut_left_od,cut_left_od_points] = find_points(exc_left_od_v,left_od_v); 
    auto imp3 = plane_generator(cut_left_od_points);
    auto [liver3,left_od_poly] = cut_run2(imp3,liver2);
    stlWriter->SetFileName("C:/code/liver_cut/new_stl/left_od.stl");
    stlWriter->SetInputData(left_od_poly);
    stlWriter->Write();
    std::cout << "Successfully saved" << std::endl;
    std::cout<<"\nleftod cut completed\n";
*/
    auto left_od_poly = liver8;
    stlWriter->SetFileName("C:/code/liver_cut/new_stl/left_od.stl");
    stlWriter->SetInputData(left_od_poly);
    stlWriter->Write();










    double mbounds[6];
    liver->GetBounds(mbounds);

    // 可视化隐式函数（等值面为0）
    vtkSmartPointer<vtkSampleFunction> sampleFunction = vtkSmartPointer<vtkSampleFunction>::New();
    sampleFunction->SetImplicitFunction(imp2); // 您的隐式函数
    sampleFunction->SetModelBounds(mbounds[0], mbounds[1], mbounds[2], mbounds[3], mbounds[4], mbounds[5]); // 设置采样范围
    sampleFunction->SetSampleDimensions(50, 50, 50); // 设置采样分辨率
    sampleFunction->ComputeNormalsOff();
    sampleFunction->Update();

    vtkSmartPointer<vtkContourFilter> contourFilter = vtkSmartPointer<vtkContourFilter>::New();
    contourFilter->SetInputConnection(sampleFunction->GetOutputPort());
    contourFilter->SetValue(0, 0.0); // 提取等值为 0 的面，即隐式函数的零等值面
    contourFilter->Update();

    vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    mapper->SetInputConnection(contourFilter->GetOutputPort());

    vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
    actor->SetMapper(mapper);
    actor->GetProperty()->SetColor(1.0, 0.0, 0.0); // 设置颜色，例如红色

    // 创建第二个映射器和演员
    vtkSmartPointer<vtkPolyDataMapper> mapper1 = vtkSmartPointer<vtkPolyDataMapper>::New();
    mapper1->SetInputData(liver);

    vtkSmartPointer<vtkActor> actor1 = vtkSmartPointer<vtkActor>::New();
    actor1->SetMapper(mapper1);
    actor1->GetProperty()->SetOpacity(0.7);  // 设置透明度
    actor1->GetProperty()->SetColor(0.0, 0.0, 1.0);

    // 创建渲染器、渲染窗口和交互器
    vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
    vtkSmartPointer<vtkRenderWindow> renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
    renderWindow->SetWindowName("STL Model and Plane Display");
    renderWindow->AddRenderer(renderer);

    vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
    renderWindowInteractor->SetRenderWindow(renderWindow);

    // 将演员添加到场景中
    renderer->AddActor(actor);
    renderer->AddActor(actor1);
    //renderer->AddActor(actor2);
    // 设置背景颜色
    renderer->SetBackground(0.1, 0.2, 0.3);

    // 调整摄像机视角，使得模型和坐标轴都能正确显示
    renderer->ResetCamera();
    renderWindow->Render();
    renderWindowInteractor->Start();


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
