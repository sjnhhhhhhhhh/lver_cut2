#ifndef CUT_RIGHT_LOBE_H
#define CUT_RIGHT_LOBE_H

#include <vtkSmartPointer.h>
#include <vtkPolyData.h>

// 声明函数来执行第二刀操作
std::tuple<vtkSmartPointer<vtkPolyData>,vtkSmartPointer<vtkPolyData>> cut_right_lobe(vtkSmartPointer<vtkPolyData> liver, 
                    vtkSmartPointer<vtkPolyData> rb1v, 
                    vtkSmartPointer<vtkPolyData> rb2v, 
                    vtkSmartPointer<vtkPolyData> rb3v, 
                    vtkSmartPointer<vtkPolyData> rf_up_v, 
                    vtkSmartPointer<vtkPolyData> rf_down_v, 
                    vtkSmartPointer<vtkPolyData> left_inside_v, 
                    vtkSmartPointer<vtkPolyData> tail);

#endif // CUT_RIGHT_LOBE_H
