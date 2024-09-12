#include <vtkSmartPointer.h>
#include <vtkSTLReader.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkPoints.h>
#include <vtkCellArray.h>
#include <vtkLine.h>
#include <vtkPolyLine.h>
#include <vtkProperty.h>
#include <fstream>
#include <sstream>
#include <iostream>

// 读取骨架多段线数据并转换为 vtkPolyData
vtkSmartPointer<vtkPolyData> ReadSkeletonFromFile(const std::string& filename) {
    std::ifstream inputFile(filename);
    if (!inputFile) {
        std::cerr << "Error: Could not open file " << filename << std::endl;
        return nullptr;
    }

    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkCellArray> lines = vtkSmartPointer<vtkCellArray>::New();
    std::string line;
    
    while (std::getline(inputFile, line)) {
        std::stringstream ss(line);
        int numPoints;
        ss >> numPoints;  // 读取多段线的点数

        vtkSmartPointer<vtkPolyLine> polyline = vtkSmartPointer<vtkPolyLine>::New();
        polyline->GetPointIds()->SetNumberOfIds(numPoints);

        for (int i = 0; i < numPoints; ++i) {
            double x, y, z;
            ss >> x >> y >> z;

            vtkIdType pointId = points->InsertNextPoint(x, y, z);
            polyline->GetPointIds()->SetId(i, pointId);
        }

        lines->InsertNextCell(polyline);
    }

    vtkSmartPointer<vtkPolyData> polyData = vtkSmartPointer<vtkPolyData>::New();
    polyData->SetPoints(points);
    polyData->SetLines(lines);

    return polyData;
}

int main() {
    // 1. 读取 STL 文件
    vtkSmartPointer<vtkSTLReader> reader = vtkSmartPointer<vtkSTLReader>::New();
    reader->SetFileName("C:/code/liver_cut/data/Hepatic_portal_vein.stl");
    reader->Update();

    vtkSmartPointer<vtkPolyDataMapper> stlMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    stlMapper->SetInputConnection(reader->GetOutputPort());

    vtkSmartPointer<vtkActor> stlActor = vtkSmartPointer<vtkActor>::New();
    stlActor->SetMapper(stlMapper);
    stlActor->GetProperty()->SetOpacity(0.5);

    // 2. 读取骨架文件
    vtkSmartPointer<vtkPolyData> skeletonPolyData = ReadSkeletonFromFile("C:/code/liver_cut/trans_data/result/skel-poly.polylines.txt");

    vtkSmartPointer<vtkPolyDataMapper> skeletonMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    skeletonMapper->SetInputData(skeletonPolyData);

    vtkSmartPointer<vtkActor> skeletonActor = vtkSmartPointer<vtkActor>::New();
    skeletonActor->SetMapper(skeletonMapper);
    skeletonActor->GetProperty()->SetColor(1.0, 0.0, 0.0);  // 红色显示骨架

    // 3. 创建渲染器和窗口
    vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
    vtkSmartPointer<vtkRenderWindow> renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
    renderWindow->AddRenderer(renderer);
    vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
    renderWindowInteractor->SetRenderWindow(renderWindow);

    // 4. 将 STL 和骨架添加到渲染器中
    renderer->AddActor(stlActor);
    renderer->AddActor(skeletonActor);
    renderer->SetBackground(0.1, 0.2, 0.4);  // 设置背景颜色

    renderWindow->Render();
    renderWindowInteractor->Start();

    return 0;
}
