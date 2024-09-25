#ifndef CUT_LEFT_INSIDE_H
#define CUT_LEFT_INSIDE_H

#include <vtkSmartPointer.h>
#include <vtkPolyData.h>

// 声明函数来执行第八刀操作
vtkSmartPointer<vtkPolyData> cut_left_inside(vtkSmartPointer<vtkPolyData> liver_left, 
                     vtkSmartPointer<vtkPolyData> left_ou_v, 
                     vtkSmartPointer<vtkPolyData> left_od_v, 
                     vtkSmartPointer<vtkPolyData> left_inside_v);

#endif // CUT_LEFT_INSIDE_H
