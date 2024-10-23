#include <iostream>
#include <unistd.h>
#include <fcntl.h>
#include <string.h>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/draw_triangulation_2.h>
#include <CGAL/point_generators_2.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Triangulation_vertex_base_2<K> Vb;
typedef CGAL::Constrained_triangulation_face_base_2<K> Fb;
typedef CGAL::Triangulation_data_structure_2<Vb,Fb> TDS;
typedef CGAL::Exact_intersections_tag Itag;
typedef CGAL::Constrained_Delaunay_triangulation_2<K, TDS, Itag> CDT;
typedef K::Point_2 Point;

#define BUFFER 1024

using namespace std;
using namespace boost::property_tree;

int fix_json() {
    string filename = "output.json";
    
    ifstream inputFile(filename);

    string fileContent((std::istreambuf_iterator<char>(inputFile)), std::istreambuf_iterator<char>());
    inputFile.close();

    std::string toReplace = "\\/";
    std::string replaceWith = "/";
    size_t pos = 0;
    
    while ((pos = fileContent.find(toReplace, pos)) != std::string::npos) {
        fileContent.replace(pos, toReplace.length(), replaceWith);
        pos += replaceWith.length();
    }

    std::ofstream outputFile(filename);

    outputFile << fileContent;
    outputFile.close();

    return 0;
}

int make_json(){
    string content_type = "CG_SHOP_2025_Solution";
    string instance_uid = "unique_instance_id";

    std::vector<string> steiner_points_x = {"123/456", "789", "1011/1213"};
    std::vector<string> steiner_points_y = {"314/159", "265", "358/979"};

    std::vector<std::pair<int, int>> edges = {{0, 7}, {7, 8}, {8, 9}};

    ptree jsonData;

    jsonData.put("content_type", content_type);

    jsonData.put("instance_uid", instance_uid);

    ptree pt_steiner_points_x;
    for (const auto& point : steiner_points_x){
        ptree temp_ptree;
        temp_ptree.put("", point);
        pt_steiner_points_x.push_back(make_pair("", temp_ptree));
    }
    jsonData.add_child("steiner_points_x", pt_steiner_points_x);

    ptree pt_steiner_points_y;
    for (const auto& point : steiner_points_y){
        ptree temp_ptree;
        temp_ptree.put("", point);
        pt_steiner_points_y.push_back(make_pair("", temp_ptree));
    }
    jsonData.add_child("steiner_points_y", pt_steiner_points_y);

    ptree pt_edges;
    for (const auto& edge : edges){
        ptree pt_edge;
        pt_edge.push_back(make_pair("", ptree(to_string(edge.first))));
        pt_edge.push_back(make_pair("", ptree(to_string(edge.second))));
        pt_edges.push_back(make_pair("", pt_edge));
    }
    jsonData.add_child("edges", pt_edges);

    try{
        write_json("output.json", jsonData);
        cout << "output.json created successfully" << endl;
    }catch (const json_parser_error &e) {
        cerr << "errot at creating JSON file: " << e.what() << endl;
        return 1;
    }

    fix_json();

    return 0;
}



CDT constrained_delaunay_function(vector<Point> points, int num_constraints, vector<vector<int>> additional_constraints){
    CDT cdt;
    for (const auto& pt : points){
        cdt.insert(pt);
    }

    for (const auto& constraint : additional_constraints){
        int x = constraint[0];
        int y = constraint[1];
        cdt.insert_constraint(points[x], points[y]);
    }

    return cdt;
}



double angle_between_vectors(const K::Vector_2& v1, const K::Vector_2& v2){
    double dot_product = v1 * v2;
    double length1 = sqrt(v1.squared_length());
    double length2 = sqrt(v2.squared_length());
    double angle = acos(dot_product / (length1 * length2)) * (180.0 / M_PI);
    return angle;
}



bool is_obtuse_triangle(const Point& p1, const Point& p2, const Point& p3){
    K::Vector_2 v1 = p2 - p1;
    K::Vector_2 v2 = p3 - p1;
    K::Vector_2 v3 = p3 - p2;

    double angle1 = angle_between_vectors(v1, v2);
    double angle2 = angle_between_vectors(-v1, v3);
    double angle3 = angle_between_vectors(v2, -v3);

    return (angle1 > 90.0 || angle2 > 90.0 || angle3 > 90.0);
}



int count_obtuse_triangles(const CDT& cdt){
    int obtuse_triangle_count = 0;
    int triangle_count = 0;

    for (auto face = cdt.finite_faces_begin(); face != cdt.finite_faces_end(); ++face){
        Point p1 = face->vertex(0)->point();
        Point p2 = face->vertex(1)->point();
        Point p3 = face->vertex(2)->point();

        triangle_count++;
        if (is_obtuse_triangle(p1, p2, p3)){
            obtuse_triangle_count++;
        }
    }

    cout << "total triangles: " << triangle_count << endl;
    return obtuse_triangle_count;
}



void flip_edge(CDT& cdt, CDT::Edge& edge) {
    CDT::Face_handle face = edge.first;
    int index = edge.second;

    Point p1 = face->vertex((index + 1) % 3)->point();
    Point p2 = face->vertex((index + 2) % 3)->point();
    Point p3 = face->vertex(index)->point();

    if (is_obtuse_triangle(p1, p2, p3)) {
        cdt.flip(face, index);
    }
}

void optimize_triangulation(CDT& cdt) {
    for (auto it = cdt.finite_edges_begin(); it != cdt.finite_edges_end(); ++it) {
        CDT::Edge edge = *it;

        CDT::Face_handle face = edge.first;
        
        Point p1 = face->vertex(0)->point();
        Point p2 = face->vertex(1)->point();
        Point p3 = face->vertex(2)->point();

        if (is_obtuse_triangle(p1, p2, p3)) {
            flip_edge(cdt, edge);
        }
    }
}




int main(void){

    char txtbuffer[BUFFER] = "";
    char filename[BUFFER] = "";

    cout << "Enter file name:" << endl;
    cin >> filename;

    ptree jsonData;
    try {
        read_json(filename, jsonData);
    } catch (const json_parser_error &e){
        cerr << "error reading file " << e.what() << endl;
        return 1;
    }

    string instance_uid = jsonData.get<string>("instance_uid");

    int num_points = jsonData.get<int>("num_points");

    vector<int> points_x;
    for (auto& item : jsonData.get_child("points_x")){
        points_x.push_back(item.second.get_value<int>());
    }

    vector<int> points_y;
    for (auto& item : jsonData.get_child("points_y")){
        points_y.push_back(item.second.get_value<int>());
    }

    vector<int> region_boundary;
    for (auto& item : jsonData.get_child("region_boundary")){
        region_boundary.push_back(item.second.get_value<int>());
    }

    int num_constraints = jsonData.get<int>("num_constraints");

    vector<vector<int>> additional_constraints;
    for (auto& item : jsonData.get_child("additional_constraints")){
        vector<int> constraint;
        for (auto& sub_item : item.second) {
            constraint.push_back(sub_item.second.get_value<int>());
        }
        additional_constraints.push_back(constraint);
    }

    vector<Point> points;
    for (int i = 0; i < num_points; ++i){
        points.push_back(Point(points_x[i], points_y[i]));
    }

    CDT cdt = constrained_delaunay_function(points, num_constraints, additional_constraints);

    //CGAL::draw(cdt);

    int obtuse_triangle_count = count_obtuse_triangles(cdt);
    cout << "Obtuse triangle count is: " << obtuse_triangle_count << endl;

    obtuse_triangle_count = count_obtuse_triangles(cdt);
    cout << "Obtuse triangle count is: " << obtuse_triangle_count << endl;

    optimize_triangulation(cdt);
    obtuse_triangle_count = count_obtuse_triangles(cdt);
    cout << "Obtuse triangle count after edge flip is: " << obtuse_triangle_count << endl;


    //make_json();
    return 0;
}