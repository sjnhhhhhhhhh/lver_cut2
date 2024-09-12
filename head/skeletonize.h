#ifndef SKELETONIZE_H
#define SKELETONIZE_H

#include <fstream>
#include <sstream>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/extract_mean_curvature_flow_skeleton.h>

typedef CGAL::Simple_cartesian<double> Kernel;
typedef Kernel::Point_3 Point;
typedef CGAL::Polyhedron_3<Kernel> Polyhedron;
typedef CGAL::Mean_curvature_flow_skeletonization<Polyhedron> Skeletonization;
typedef Skeletonization::Skeleton Skeleton;

// 加载并检查模型是否是三角形网格
bool load_and_check_model(const std::string& filepath, Polyhedron& mesh);

// 提取模型的骨架
void extract_skeleton(const Polyhedron& mesh, Skeleton& skeleton);

// 用于输出骨架的多段线
struct Display_polylines {
    const Skeleton& skeleton;
    std::ofstream& out;
    int polyline_size;
    std::stringstream sstr;
    Display_polylines(const Skeleton& skeleton, std::ofstream& out);
    void start_new_polyline();
    void add_node(Skeleton::vertex_descriptor v);
    void end_polyline();
};

// 输出骨架的多段线到文件
void output_skeleton_polylines(const Skeleton& skeleton, const std::string& filepath);

// 输出骨架点和对应的曲面点到文件
void output_skeleton_correspondence(const Skeleton& skeleton, const Polyhedron& mesh, const std::string& filepath);




#endif
