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
typedef CDT::Face_handle Face_handle;

#define BUFFER 1024

using namespace std;
using namespace boost::property_tree;

int replace_json(string oldText, string newText) {
    int result = 1;
    string filename = "output.json";
    
    ifstream inputFile(filename);

    string fileContent((istreambuf_iterator<char>(inputFile)), istreambuf_iterator<char>());
    inputFile.close();

    size_t pos = 0;
    
    while ((pos = fileContent.find(oldText, pos)) != string::npos) {
        fileContent.replace(pos, oldText.length(), newText);
        pos += newText.length();
        result = 0;
    }

    std::ofstream outputFile(filename);

    outputFile << fileContent;
    outputFile.close();

    return result;
}

int make_json(vector<Point> steiner_points){
    string content_type = "CG_SHOP_2025_Solution";
    string instance_uid = "unique_instance_id";

    vector<pair<int, int>> edges = {{0, 7}, {7, 8}, {8, 9}};

    ptree jsonData;

    jsonData.put("content_type", content_type);

    jsonData.put("instance_uid", instance_uid);

    ptree pt_steiner_points_x;
    for (const auto& point : steiner_points){
        ptree temp_ptree;
        temp_ptree.put("", point[0]);
        pt_steiner_points_x.push_back(make_pair("", temp_ptree));
    }
    jsonData.add_child("steiner_points_x", pt_steiner_points_x);

    ptree pt_steiner_points_y;
    for (const auto& point : steiner_points){
        ptree temp_ptree;
        temp_ptree.put("", point[1]);
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

    //replace_json("\\/", "/");

    return 0;
}



CDT constrained_delaunay_function(vector<Point> points, int num_constraints, vector<vector<int>> additional_constraints){
    CDT cdt;
    for (const auto& pt : points){
        cdt.insert(pt);
    }

    for (const auto& constraint : additional_constraints){
        int a = constraint[0];
        int b = constraint[1];
        cdt.insert_constraint(points[a], points[b]);
    }

    return cdt;
}



bool is_obtuse_triangle(const Point& p1, const Point& p2, const Point& p3){
    K::Vector_2 v1 = p2 - p1;
    K::Vector_2 v2 = p3 - p1;
    K::Vector_2 v3 = p3 - p2;

    double dot1 = v1 * v2;
    double dot2 = -v1 * v3;
    double dot3 = v2 * -v3;

    return (dot1 < 0 || dot2 < 0 || dot3 < 0);
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



pair<Point, Point> find_longest_edge(const CDT& cdt){
    pair<Point, Point> longest_edge;
    double max_distance = 0;

    for (auto edge = cdt.finite_edges_begin(); edge != cdt.finite_edges_end(); ++edge){
        Point p1 = cdt.segment(*edge).source();
        Point p2 = cdt.segment(*edge).target();

        double distance = CGAL::squared_distance(p1, p2);

        if (distance > max_distance) {
            max_distance = distance;
            longest_edge = make_pair(p1, p2);
        }
    }

    return longest_edge;
}



Face_handle find_obtuse_triangle(CDT& cdt){
    for (auto face = cdt.finite_faces_begin(); face != cdt.finite_faces_end(); ++face) {
        Point p1 = face->vertex(0)->point();
        Point p2 = face->vertex(1)->point();
        Point p3 = face->vertex(2)->point();

        if (is_obtuse_triangle(p1, p2, p3)) {
            return face;
        }
    }

    return nullptr;
}



Point steiner_at_midpoint (CDT& cdt){
    pair<Point, Point> edge = find_longest_edge(cdt);

    Point p1 = edge.first;
    Point p2 = edge.second;

    Point steiner_point = CGAL::midpoint(p1,p2);
    cdt.insert(steiner_point);
    return steiner_point;
}




Point steiner_at_circumcenter(CDT& cdt){
    Face_handle triangle = find_obtuse_triangle(cdt);
    if (triangle == nullptr){
        return Point(nan(""), nan(""));
    }

    Point p1 = triangle->vertex(0)->point();
    Point p2 = triangle->vertex(1)->point();
    Point p3 = triangle->vertex(2)->point();

    Point steiner_point = CGAL::circumcenter(p1, p2, p3);

    cdt.insert(steiner_point);
    return steiner_point;
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

    vector<Point> points;
    for (int i = 0; i < num_points; ++i){
        points.push_back(Point(points_x[i], points_y[i]));
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

    CDT cdt = constrained_delaunay_function(points, num_constraints, additional_constraints);

    CGAL::draw(cdt);
    int obtuse_triangle_count = count_obtuse_triangles(cdt);
    cout << "Obtuse triangle count is: " << obtuse_triangle_count << endl;

    vector<Point> steiner_points;

    int method;
    cout << "Enter method number (1,2,3)" << endl;
    cin >> method;

    for (int i = 0; i <10; i++){
        Point steiner_point;
        if (method == 1) steiner_point = steiner_at_midpoint(cdt);
        else if (method == 2) steiner_point = steiner_at_circumcenter(cdt);
        else if (method == 3); // Point steiner_point = steiner_at_3(cdt);
        else {
            cout << "Wrong method imput." << endl;
            return 2;
        }
        if (steiner_point[0] != nan("") && steiner_point[1] != nan("")) steiner_points.push_back(steiner_point);
    }
    CGAL::draw(cdt);
    obtuse_triangle_count = count_obtuse_triangles(cdt);
    cout << "Obtuse triangle count is: " << obtuse_triangle_count << endl;

    make_json(steiner_points);

    return 0;
}