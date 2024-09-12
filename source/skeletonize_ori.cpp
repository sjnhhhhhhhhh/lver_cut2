#include <list>
#include <cassert>



#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/extract_mean_curvature_flow_skeleton.h>
#include <CGAL/boost/graph/split_graph_into_polylines.h>


typedef CGAL::Simple_cartesian<double>                        Kernel;
typedef Kernel::Point_3                                       Point;
typedef CGAL::Polyhedron_3<Kernel>                            Polyhedron;
typedef boost::graph_traits<Polyhedron>::vertex_descriptor    vertex_descriptor;
typedef CGAL::Mean_curvature_flow_skeletonization<Polyhedron> Skeletonization;
typedef Skeletonization::Skeleton                             Skeleton;
typedef Skeleton::vertex_descriptor                           Skeleton_vertex;
typedef Skeleton::edge_descriptor                             Skeleton_edge;


//仅需要最大多段线的显示
struct Display_polylines {
    const Skeleton& skeleton;
    std::ofstream& out;
    int polyline_size;
    std::stringstream sstr;
    Display_polylines(const Skeleton& skeleton, std::ofstream& out)
        : skeleton(skeleton), out(out)
    {}
    void start_new_polyline() {
        polyline_size = 0;
        sstr.str("");
        sstr.clear();
    }
    void add_node(Skeleton_vertex v) {
        ++polyline_size;
        sstr << " " << skeleton[v].point;
    }
    void end_polyline()
    {
        out << polyline_size << sstr.str() << "\n";
    }
};

int main(int argc, const char** argv)
{
    //****************************************获取数据*****************************************************
    std::ifstream input((argc > 1) ? argv[1] : CGAL::data_file_path("C:/code/liver_cut/trans_data/Hepatic_portal_vein.off"));
    Polyhedron tmesh;
    input >> tmesh;
    if (!CGAL::is_triangle_mesh(tmesh))
    {
        std::cout << "Input geometry is not triangulated." << std::endl;
        return EXIT_FAILURE;
    }

    //***************************************网格骨架提取*************************************************
    Skeleton skeleton;
    CGAL::extract_mean_curvature_flow_skeleton(tmesh, skeleton);
    std::cout << "Number of vertices of the skeleton: " << boost::num_vertices(skeleton) << "\n";
    std::cout << "Number of edges of the skeleton: " << boost::num_edges(skeleton) << "\n";

    //****************************************输出结果*************************************************
    //输出骨架所有的边
    std::ofstream output("C:/code/liver_cut/trans_data/result/skel-poly.polylines.txt");
    Display_polylines display(skeleton, output);
    CGAL::split_graph_into_polylines(skeleton, display);
    output.close();

    //输出骨架所有的点以及其相对应的曲面点
    output.open("C:/code/liver_cut/trans_data/result/correspondance-poly.polylines.txt");
    for (Skeleton_vertex v : CGAL::make_range(vertices(skeleton)))
        for (vertex_descriptor vd : skeleton[v].vertices)
            output << /*"2 " <<*/ skeleton[v].point << " "
            << get(CGAL::vertex_point, tmesh, vd) << "\n";
    output.close();
    //std::cout <<"计算结束" << std::endl;
    //system("pause");




    return EXIT_SUCCESS;
}

