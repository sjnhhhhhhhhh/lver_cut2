# CMakeLists.txt: test1 的 CMake 项目，在此处包括源代码并定义项目特定的逻辑。

# 设置源文件列表，包括主程序和模块
set(SOURCE_FILES
    main.cpp       # 主程序文件，包含main函数
    source/liver_cut.cpp   # keypoint模块实现文件
    source/cut_tail.cpp
    source/cut_right_lobe.cpp
    source/cut_left_lobe.cpp
    source/cut_rb2b.cpp
    source/cut_rb3b.cpp
    source/cut_rf_up.cpp
    source/cut_rb1b_rf_down.cpp
    source/cut_left_inside.cpp
    source/cut_left_ou_od.cpp
    #keypointfinder.cpp  # 其他模块实现文件
    #keypointfinder2.cpp
    #source/STL2OFF.cpp    # 其他模块实现文件
    #source/skeletonize.cpp # 其他模块实现文件
    #test.cpp
    #testopen3d.cpp
    #cut_test.cpp
)

# 生成可执行文件，并使用源文件列表
add_executable(liver_cut ${SOURCE_FILES})


if (CMAKE_VERSION VERSION_GREATER 3.12)
  set_property(TARGET liver_cut PROPERTY CXX_STANDARD 20)
endif()

include_directories("C:/source/InsightToolkit-5.1.0/Modules/Bridge/VtkGlue/include")
include_directories(${CMAKE_SOURCE_DIR}/src/head)

# 手动添加 pybind11 头文件路径
include_directories("C:/Program Files/Python311/Lib/site-packages/pybind11/include")
# 链接 Python 解释器
find_package(Python3 COMPONENTS Interpreter Development REQUIRED)
target_link_libraries(liver_cut PRIVATE Python3::Python)


# 查找并链接VTK库
find_package(VTK REQUIRED)
include(${VTK_USE_FILE}) # 确保 VTK 配置正确
target_link_libraries(liver_cut PRIVATE ${VTK_LIBRARIES})

# 查找并链接ITK库
set(CMAKE_PREFIX_PATH "C:/source/ITK" ${CMAKE_PREFIX_PATH})
find_package(ITK REQUIRED)
include(${ITK_USE_FILE})
target_link_libraries(liver_cut PRIVATE ${ITK_LIBRARIES})

# 查找CGAL库
find_package(CGAL CONFIG REQUIRED)
# 添加可执行文件并链接CGAL库
target_link_libraries(liver_cut PRIVATE CGAL::CGAL)

include("C:/source/eigen-3.4.0/cmake/UseEigen3.cmake")

# 查找PCL库
find_package(PCL CONFIG REQUIRED)
target_link_libraries(liver_cut PRIVATE ${PCL_LIBRARIES})

# 添加Open3D根目录路径到CMake的搜索路径
list(APPEND CMAKE_PREFIX_PATH "C:/source/Open3d/open3d-devel-windows-amd64-0.17.0")

find_package(Open3D REQUIRED)

# 链接Open3D库
target_link_libraries(liver_cut PRIVATE Open3D::Open3D)

# 由于下载的Open3D是动态库，需要将Open3D.dll拷贝到执行文件的目录下，以下自定义指令完成了这一工作。
add_custom_command(TARGET liver_cut POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy
    "${Open3D_DIR}/../bin/Open3D.dll"
    $<TARGET_FILE_DIR:liver_cut>)



# 在目标创建后使用 vtk_module_autoinit
if (VTK_VERSION VERSION_GREATER "9.00.0")
    vtk_module_autoinit(
        TARGETS liver_cut
        MODULES ${VTK_LIBRARIES}
    )
endif ()

include_directories(${PROJECT_SOURCE_DIR}/include)
include_directories(${EIGEN3_INCLUDE_DIR})