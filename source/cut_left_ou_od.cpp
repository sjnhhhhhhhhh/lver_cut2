#include "cut_left_ou_od.h"
#include "liver_cut.h"  // 包含你已经模块化的函数，比如 plane_generator, cut_run2 等
#include <vtkSTLWriter.h>
#include <iostream>

// 定义函数，执行第九刀的操作
void cut_left_ou_od(vtkSmartPointer<vtkPolyData> liver_l1, 
                    vtkSmartPointer<vtkPolyData> left_od_v, 
                    vtkSmartPointer<vtkPolyData> left_ou_v) 
{
    // 第9刀，左后上段和下段
    auto exc_left_ou_v = mergePolyData(left_od_v);
    auto [cut_left_ou, cut_left_ou_points] = find_points(exc_left_ou_v, left_ou_v);
    auto imp9 = plane_generator(cut_left_ou);
    auto [left_od_poly, left_ou_poly] = cut_run2(imp9, liver_l1);

    // 保存左后上段切割后的 STL 文件
    vtkSmartPointer<vtkSTLWriter> stlWriter = vtkSmartPointer<vtkSTLWriter>::New();
    stlWriter->SetFileName("C:/code/liver_cut/new_stl/left_ou.stl");
    stlWriter->SetInputData(left_ou_poly);
    stlWriter->Write();
    std::cout << "left_ou cut completed" << std::endl;

    // 保存左后下段切割后的 STL 文件
    stlWriter->SetFileName("C:/code/liver_cut/new_stl/left_od.stl");
    stlWriter->SetInputData(left_od_poly);
    stlWriter->Write();
    std::cout << "left_od cut completed" << std::endl;
}
