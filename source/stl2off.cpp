#include "STL2OFF.h"
#include <vtkSmartPointer.h>
#include <vtkSTLReader.h>
#include <vtkPolyData.h>
#include <vtkCell.h>
#include <fstream>
#include <iostream>

void STL2OFF(const std::string& stl_filename, const std::string& off_filename) {
    // 读取 STL 文件
    vtkSmartPointer<vtkSTLReader> reader = vtkSmartPointer<vtkSTLReader>::New();
    reader->SetFileName(stl_filename.c_str());
    reader->Update();

    vtkSmartPointer<vtkPolyData> polyData = reader->GetOutput();

    // 打开输出文件
    std::ofstream off_file(off_filename);
    if (!off_file.is_open()) {
        std::cerr << "Error: Unable to open OFF file for writing." << std::endl;
        return;
    }

    // 写入 OFF 文件头
    off_file << "OFF\n";
    off_file << polyData->GetNumberOfPoints() << " " << polyData->GetNumberOfCells() << " 0\n";

    // 写入顶点数据
    double point[3];
    for (vtkIdType i = 0; i < polyData->GetNumberOfPoints(); i++) {
        polyData->GetPoint(i, point);
        off_file << point[0] << " " << point[1] << " " << point[2] << "\n";
    }

    // 写入面数据
    for (vtkIdType i = 0; i < polyData->GetNumberOfCells(); i++) {
        vtkCell* cell = polyData->GetCell(i);
        if (cell->GetNumberOfPoints() == 3) { // 确保是三角形
            off_file << "3 "; // 每个面都有三个顶点
            for (vtkIdType j = 0; j < cell->GetNumberOfPoints(); j++) {
                off_file << cell->GetPointId(j) << " ";
            }
            off_file << "\n";
        }
    }

    off_file.close();
    std::cout << "Conversion to OFF completed successfully!" << std::endl;
}
