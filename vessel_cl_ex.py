import vtk
import sys
# 读取分割的血管图像
reader = vtk.vtkXMLPolyDataReader()
reader.SetFileName('Hepatic_portal_vein.stl')
reader.Update()

# 提取中心线
centerline_filter = vmtkscripts.vmtkCenterlines()
centerline_filter.Surface = reader.GetOutput()
centerline_filter.Execute()

# 保存中心线到文件
writer = vtk.vtkXMLPolyDataWriter()
writer.SetFileName('centerlines.vtp')
writer.SetInputData(centerline_filter.Centerlines)
writer.Write()
