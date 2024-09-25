#include "cut_tail.h"
#include "liver_cut.h" 
#include <vtkSTLWriter.h>
#include <iostream>

// 定义函数，执行第一刀的操作
vtkSmartPointer<vtkPolyData> cut_tail(vtkSmartPointer<vtkPolyData> liver, 
              vtkSmartPointer<vtkPolyData> rb1v, 
              vtkSmartPointer<vtkPolyData> rb2v, 
              vtkSmartPointer<vtkPolyData> rb3v, 
              vtkSmartPointer<vtkPolyData> rf_up_v, 
              vtkSmartPointer<vtkPolyData> rf_down_v, 
              vtkSmartPointer<vtkPolyData> left_inside_v, 
              vtkSmartPointer<vtkPolyData> left_ou_v, 
              vtkSmartPointer<vtkPolyData> left_od_v, 
              vtkSmartPointer<vtkPolyData> tail) 
{
    // 第1刀，尾状段
    auto exc_tail = mergePolyData(rb1v, rb2v, rb3v, rf_up_v, rf_down_v, left_inside_v, left_ou_v, left_od_v);
    auto [cut_tail, cut_tail_points] = find_points(exc_tail, tail);
    auto imp1 = plane_generator(cut_tail);
    auto [liver1, tail_poly] = cut_run2(imp1, liver);

    // 保存切割后的 STL 文件
    vtkSmartPointer<vtkSTLWriter> stlWriter = vtkSmartPointer<vtkSTLWriter>::New();
    stlWriter->SetFileName("C:/code/liver_cut/new_stl/tail.stl");
    stlWriter->SetInputData(tail_poly);
    stlWriter->Write();

    std::cout << "Successfully saved tail cut" << std::endl;
    std::cout << "Tail cut completed" << std::endl;
    return liver1;
}
