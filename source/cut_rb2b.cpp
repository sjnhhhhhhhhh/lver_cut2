#include "cut_rb2b.h"
#include "liver_cut.h"  // 包含你已经模块化的函数，比如 plane_generator, cut_run2 等
#include <vtkSTLWriter.h>
#include <iostream>

// 定义函数，执行第五刀的操作
vtkSmartPointer<vtkPolyData> cut_rb2b(vtkSmartPointer<vtkPolyData> liver_r1, 
              vtkSmartPointer<vtkPolyData> rb1v, 
              vtkSmartPointer<vtkPolyData> rb3v, 
              vtkSmartPointer<vtkPolyData> rf_up_v, 
              vtkSmartPointer<vtkPolyData> tail, 
              vtkSmartPointer<vtkPolyData> rb2v) 
{
    // 第5刀，右后2段
    auto exc_rb2v = mergePolyData(rb1v, rb3v, rf_up_v, tail);
    auto [rb2b, rb2b_points] = find_points(exc_rb2v, rb2v);
    auto imp5 = plane_generator(rb2b);
    auto [liver_r2, rb2b_poly] = cut_run2(imp5, liver_r1);

    // 保存右后2段切割后的 STL 文件
    vtkSmartPointer<vtkSTLWriter> stlWriter = vtkSmartPointer<vtkSTLWriter>::New();
    stlWriter->SetFileName("C:/code/liver_cut/new_stl/rb2b.stl");
    stlWriter->SetInputData(rb2b_poly);
    stlWriter->Write();

    std::cout << "Successfully saved rb2b cut" << std::endl;
    std::cout << "rb2b cut completed" << std::endl;
    
    return liver_r2;
}
