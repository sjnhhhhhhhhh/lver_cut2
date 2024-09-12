#include "skeletonize.h"
#include <iostream>
#include <fstream>
#include <CGAL/boost/graph/split_graph_into_polylines.h>
//****************************************获取数据*****************************************************
bool load_and_check_model(const std::string& filepath, Polyhedron& mesh) {
    std::cout<<"obtaining data\n";
    std::ifstream input(filepath);
    input >> mesh;
    if (!CGAL::is_triangle_mesh(mesh)) {
        std::cerr << "Input geometry is not triangulated." << std::endl;
        return false;
        std::cout<<"obtaining fail\n";
    }
    return true;
    std::cout<<"obtaining success\n";
}
//***************************************网格骨架提取*************************************************
void extract_skeleton(const Polyhedron& mesh, Skeleton& skeleton) {
    std::cout<<"extracting skeleton\n";
    CGAL::extract_mean_curvature_flow_skeleton(mesh, skeleton);
    std::cout << "Number of vertices of the skeleton: " << boost::num_vertices(skeleton) << "\n";
    std::cout << "Number of edges of the skeleton: " << boost::num_edges(skeleton) << "\n";
}



Display_polylines::Display_polylines(const Skeleton& skeleton, std::ofstream& out)
    : skeleton(skeleton), out(out) {}

void Display_polylines::start_new_polyline() {
    polyline_size = 0;
    sstr.str("");
    sstr.clear();
}

void Display_polylines::add_node(Skeleton::vertex_descriptor v) {
    ++polyline_size;
    sstr << " " << skeleton[v].point;
}

void Display_polylines::end_polyline() {
    out << polyline_size << sstr.str() << "\n";
}

void output_skeleton_polylines(const Skeleton& skeleton, const std::string& filepath) {
    std::ofstream output(filepath);
    Display_polylines display(skeleton, output);
    CGAL::split_graph_into_polylines(skeleton, display);
    output.close();
}

void output_skeleton_correspondence(const Skeleton& skeleton, const Polyhedron& mesh, const std::string& filepath) {
    std::ofstream output(filepath);
    for (auto v : CGAL::make_range(vertices(skeleton)))
        for (auto vd : skeleton[v].vertices)
            output << skeleton[v].point << " " << get(CGAL::vertex_point, mesh, vd) << "\n";
    output.close();
}