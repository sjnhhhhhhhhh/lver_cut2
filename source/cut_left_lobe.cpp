#include "cut_left_lobe.h"
#include "liver_cut.h"  // 包含你已经模块化的函数，比如 plane_generator, cut_run2 等
#include <vtkSTLWriter.h>
#include <iostream>

// 定义函数，执行第三刀的操作
vtkSmartPointer<vtkPolyData> cut_left_lobe(vtkSmartPointer<vtkPolyData> liver_left_pre, 
                   vtkSmartPointer<vtkPolyData> left_inside_v, 
                   vtkSmartPointer<vtkPolyData> left_od_v, 
                   vtkSmartPointer<vtkPolyData> left_ou_v, 
                   vtkSmartPointer<vtkPolyData> tail) 
{
    // 第3刀，提取左半叶
    auto left_side_v = mergePolyData(left_inside_v, left_od_v, left_ou_v);
    auto [cut_left_side, cut_left_side_points] = find_points(left_side_v, tail);
    auto imp3 = plane_generator(cut_left_side);
    auto [liver_left, t] = cut_run2(imp3, liver_left_pre);

    // 保存左半叶切割后的 STL 文件
    vtkSmartPointer<vtkSTLWriter> stlWriter = vtkSmartPointer<vtkSTLWriter>::New();
    stlWriter->SetFileName("C:/code/liver_cut/new_stl/liver_left.stl");
    stlWriter->SetInputData(liver_left);
    stlWriter->Write();

    std::cout << "Successfully saved liver_left cut" << std::endl;
    std::cout << "Liver_left cut completed" << std::endl;

    return liver_left;
}
