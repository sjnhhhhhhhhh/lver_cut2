#ifndef CUT_LEFT_OU_OD_H
#define CUT_LEFT_OU_OD_H

#include <vtkSmartPointer.h>
#include <vtkPolyData.h>

// 声明函数来执行第九刀操作
void cut_left_ou_od(vtkSmartPointer<vtkPolyData> liver_l1, 
                    vtkSmartPointer<vtkPolyData> left_od_v, 
                    vtkSmartPointer<vtkPolyData> left_ou_v);

#endif // CUT_LEFT_OU_OD_H
