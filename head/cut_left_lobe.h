#ifndef CUT_LEFT_LOBE_H
#define CUT_LEFT_LOBE_H

#include <vtkSmartPointer.h>
#include <vtkPolyData.h>

// 声明函数来执行第三刀操作
vtkSmartPointer<vtkPolyData> cut_left_lobe(vtkSmartPointer<vtkPolyData> liver_left_pre, 
                   vtkSmartPointer<vtkPolyData> left_inside_v, 
                   vtkSmartPointer<vtkPolyData> left_od_v, 
                   vtkSmartPointer<vtkPolyData> left_ou_v, 
                   vtkSmartPointer<vtkPolyData> tail);

#endif // CUT_LEFT_LOBE_H
