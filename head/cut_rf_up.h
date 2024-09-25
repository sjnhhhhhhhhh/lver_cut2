#ifndef CUT_RF_UP_H
#define CUT_RF_UP_H

#include <vtkSmartPointer.h>
#include <vtkPolyData.h>

// 声明函数来执行第六刀操作
vtkSmartPointer<vtkPolyData> cut_rf_up(vtkSmartPointer<vtkPolyData> liver_r2, 
               vtkSmartPointer<vtkPolyData> rb1v, 
               vtkSmartPointer<vtkPolyData> rf_down_v, 
               vtkSmartPointer<vtkPolyData> rf_up_v);

#endif // CUT_RF_UP_H
