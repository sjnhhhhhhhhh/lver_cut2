#ifndef CUT_TAIL_H
#define CUT_TAIL_H

#include <vtkSmartPointer.h>
#include <vtkPolyData.h>

// 声明函数来执行第一刀操作
vtkSmartPointer<vtkPolyData> cut_tail(vtkSmartPointer<vtkPolyData> liver, 
              vtkSmartPointer<vtkPolyData> rb1v, 
              vtkSmartPointer<vtkPolyData> rb2v, 
              vtkSmartPointer<vtkPolyData> rb3v, 
              vtkSmartPointer<vtkPolyData> rf_up_v, 
              vtkSmartPointer<vtkPolyData> rf_down_v, 
              vtkSmartPointer<vtkPolyData> left_inside_v, 
              vtkSmartPointer<vtkPolyData> left_ou_v, 
              vtkSmartPointer<vtkPolyData> left_od_v, 
              vtkSmartPointer<vtkPolyData> tail);

#endif // CUT_TAIL_H
