#include "cut_rb3b.h"
#include "liver_cut.h" 
#include <vtkSTLWriter.h>
#include <iostream>

// 定义函数，执行第四刀的操作
vtkSmartPointer<vtkPolyData> cut_rb3b(vtkSmartPointer<vtkPolyData> liver_right, 
              vtkSmartPointer<vtkPolyData> rb1v, 
              vtkSmartPointer<vtkPolyData> rb2v, 
              vtkSmartPointer<vtkPolyData> rf_up_v, 
              vtkSmartPointer<vtkPolyData> rf_down_v, 
              vtkSmartPointer<vtkPolyData> rb3v, 
              vtkSmartPointer<vtkPolyData> tail) 
{
    // 第4刀，右后3段
    auto exc_rb3v = mergePolyData(rb1v, rb2v, rf_down_v, rf_up_v, tail);
    auto [rb3b, rb3b_points] = find_points(exc_rb3v, rb3v);
    auto imp4 = plane_generator(rb3b);
    auto [liver_r1, rb3b_poly] = cut_run2(imp4, liver_right);

    // 保存右后3段切割后的 STL 文件
    vtkSmartPointer<vtkSTLWriter> stlWriter = vtkSmartPointer<vtkSTLWriter>::New();
    stlWriter->SetFileName("C:/code/liver_cut/new_stl/rb3b.stl");
    stlWriter->SetInputData(rb3b_poly);
    stlWriter->Write();

    std::cout << "Successfully saved rb3b cut" << std::endl;
    std::cout << "rb3b cut completed" << std::endl;

    return liver_r1;
}
