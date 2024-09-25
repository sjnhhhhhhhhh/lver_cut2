#ifndef CUT_RB2B_H
#define CUT_RB2B_H

#include <vtkSmartPointer.h>
#include <vtkPolyData.h>

// 声明函数来执行第五刀操作
vtkSmartPointer<vtkPolyData> cut_rb2b(vtkSmartPointer<vtkPolyData> liver_r1, 
              vtkSmartPointer<vtkPolyData> rb1v, 
              vtkSmartPointer<vtkPolyData> rb3v, 
              vtkSmartPointer<vtkPolyData> rf_up_v, 
              vtkSmartPointer<vtkPolyData> tail, 
              vtkSmartPointer<vtkPolyData> rb2v);

#endif // CUT_RB2B_H
