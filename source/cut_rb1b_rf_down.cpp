#include "cut_rb1b_rf_down.h"
#include "liver_cut.h"  // 包含你已经模块化的函数，比如 plane_generator, cut_run2 等
#include <vtkSTLWriter.h>
#include <iostream>

// 定义函数，执行第七刀的操作
void cut_rb1b_rf_down(vtkSmartPointer<vtkPolyData> liver_r3, 
                      vtkSmartPointer<vtkPolyData> rb1v, 
                      vtkSmartPointer<vtkPolyData> rf_down_v) 
{
    // 第7刀，右后1段与右前下段
    auto [rb1b, rb1b_points] = find_points(rf_down_v, rb1v);
    auto imp7 = plane_generator(rb1b);
    auto [rf_down_poly, rb1b_poly] = cut_run2(imp7, liver_r3);

    // 保存右后1段切割后的 STL 文件
    vtkSmartPointer<vtkSTLWriter> stlWriter = vtkSmartPointer<vtkSTLWriter>::New();
    stlWriter->SetFileName("C:/code/liver_cut/new_stl/rb1b.stl");
    stlWriter->SetInputData(rb1b_poly);
    stlWriter->Write();
    std::cout << "rb1b cut completed" << std::endl;

    // 保存右前下段切割后的 STL 文件
    stlWriter->SetFileName("C:/code/liver_cut/new_stl/rf_down.stl");
    stlWriter->SetInputData(rf_down_poly);
    stlWriter->Write();
    std::cout << "rf_down cut completed" << std::endl;

}
