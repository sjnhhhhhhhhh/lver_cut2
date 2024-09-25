#include "cut_rf_up.h"
#include "liver_cut.h"  // 包含你已经模块化的函数，比如 plane_generator, cut_run2 等
#include <vtkSTLWriter.h>
#include <iostream>

// 定义函数，执行第六刀的操作
vtkSmartPointer<vtkPolyData> cut_rf_up(vtkSmartPointer<vtkPolyData> liver_r2, 
               vtkSmartPointer<vtkPolyData> rb1v, 
               vtkSmartPointer<vtkPolyData> rf_down_v, 
               vtkSmartPointer<vtkPolyData> rf_up_v) 
{
    // 第6刀，右前上段
    auto exc_rf_up_v = mergePolyData(rf_down_v, rb1v);
    auto [rf_up_b, rf_up_v_points] = find_points(exc_rf_up_v, rf_up_v);
    auto imp6 = plane_generator(rf_up_b);
    auto [liver_r3, rf_up_v_poly] = cut_run2(imp6, liver_r2);

    // 保存右前上段切割后的 STL 文件
    vtkSmartPointer<vtkSTLWriter> stlWriter = vtkSmartPointer<vtkSTLWriter>::New();
    stlWriter->SetFileName("C:/code/liver_cut/new_stl/rf_up.stl");
    stlWriter->SetInputData(rf_up_v_poly);
    stlWriter->Write();

    std::cout << "Successfully saved rf_up cut" << std::endl;
    std::cout << "rf_up cut completed" << std::endl;

    return liver_r3;
}
