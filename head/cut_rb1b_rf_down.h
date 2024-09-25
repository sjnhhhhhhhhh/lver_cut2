#ifndef CUT_RB1B_RF_DOWN_H
#define CUT_RB1B_RF_DOWN_H

#include <vtkSmartPointer.h>
#include <vtkPolyData.h>

// 声明函数来执行第七刀操作
void cut_rb1b_rf_down(vtkSmartPointer<vtkPolyData> liver_r3, 
                      vtkSmartPointer<vtkPolyData> rb1v, 
                      vtkSmartPointer<vtkPolyData> rf_down_v);

#endif // CUT_RB1B_RF_DOWN_H
