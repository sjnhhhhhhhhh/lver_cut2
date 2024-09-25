#include "cut_left_inside.h"
#include "liver_cut.h"  // 包含你已经模块化的函数，比如 plane_generator, cut_run2 等
#include <vtkSTLWriter.h>
#include <iostream>

// 定义函数，执行第八刀的操作
vtkSmartPointer<vtkPolyData> cut_left_inside(vtkSmartPointer<vtkPolyData> liver_left, 
                     vtkSmartPointer<vtkPolyData> left_ou_v, 
                     vtkSmartPointer<vtkPolyData> left_od_v, 
                     vtkSmartPointer<vtkPolyData> left_inside_v) 
{
    // 第8刀，左内叶段
    auto exc_left_inside_v1 = mergePolyData(left_ou_v, left_od_v);
    auto [cut_left_inside1, cut_left_inside_points1] = find_points(exc_left_inside_v1, left_inside_v);
    auto imp8 = plane_generator(cut_left_inside1);
    auto [left_inside_poly, liver_l1] = cut_run2(imp8, liver_left);

    // 保存左内叶段切割后的 STL 文件
    vtkSmartPointer<vtkSTLWriter> stlWriter = vtkSmartPointer<vtkSTLWriter>::New();
    stlWriter->SetFileName("C:/code/liver_cut/new_stl/left_inside.stl");
    stlWriter->SetInputData(left_inside_poly);
    stlWriter->Write();

    std::cout << "Successfully saved left_inside cut" << std::endl;
    std::cout << "left_inside cut completed" << std::endl;

    return liver_l1;
}
