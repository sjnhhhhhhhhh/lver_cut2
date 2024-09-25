#ifndef CUT_RB3B_H
#define CUT_RB3B_H

#include <vtkSmartPointer.h>
#include <vtkPolyData.h>

// 声明函数来执行第四刀操作
vtkSmartPointer<vtkPolyData> cut_rb3b(vtkSmartPointer<vtkPolyData> liver_right, 
              vtkSmartPointer<vtkPolyData> rb1v, 
              vtkSmartPointer<vtkPolyData> rb2v, 
              vtkSmartPointer<vtkPolyData> rf_up_v, 
              vtkSmartPointer<vtkPolyData> rf_down_v, 
              vtkSmartPointer<vtkPolyData> rb3v, 
              vtkSmartPointer<vtkPolyData> tail);

#endif // CUT_RB3B_H
