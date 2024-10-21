#include <iostream>
#include <unistd.h>
#include <fcntl.h>
#include <string.h>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
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
    for (const auto& point : steiner_points_x) {
        ptree temp_ptree;
        temp_ptree.put("", point);
        pt_steiner_points_x.push_back(make_pair("", temp_ptree));
    }
    jsonData.add_child("steiner_points_x", pt_steiner_points_x);

    ptree pt_steiner_points_y;
    for (const auto& point : steiner_points_y) {
        ptree temp_ptree;
        temp_ptree.put("", point);
        pt_steiner_points_y.push_back(make_pair("", temp_ptree));
    }
    jsonData.add_child("steiner_points_y", pt_steiner_points_y);

    ptree pt_edges;
    for (const auto& edge : edges) {
        ptree pt_edge;
        pt_edge.push_back(make_pair("", ptree(to_string(edge.first))));
        pt_edge.push_back(make_pair("", ptree(to_string(edge.second))));
        pt_edges.push_back(make_pair("", pt_edge));
    }
    jsonData.add_child("edges", pt_edges);

    try {
        write_json("output.json", jsonData);
        cout << "output.json created successfully" << endl;
    } catch (const json_parser_error &e) {
        cerr << "errot at creating JSON file: " << e.what() << endl;
        return 1;
    }

    return 0;
}

int constrained_delaunay_function(vector<Point> points, int num_constraints, vector<vector<int>> additional_constraints){
    CDT cdt;
    for (const auto& pt : points) {
        cdt.insert(pt);
    }

    for (int i = 0; i<num_constraints ; i++){
        cdt.insert_constraint(additional_constraints[i][0], additional_constraints[i][1]);
    }

    for (auto edge = cdt.edges_begin(); edge != cdt.edges_end(); ++edge) {
        CDT::Segment segment = cdt.segment(*edge);
        std::cout << segment << std::endl;
    }

    return 0;
}


int main(void){

    char txtbuffer[BUFFER] = "";
    char filename[BUFFER] = "";

    cout << "Enter file name:" << endl;
    cin >> filename;

    ptree jsonData;
    try {
        read_json(filename, jsonData);
    } catch (const json_parser_error &e) {
        cerr << "error reading file " << e.what() << endl;
        return 1;
    }

    string instance_uid = jsonData.get<string>("instance_uid");

    int num_points = jsonData.get<int>("num_points");

    vector<int> points_x;
    for (auto& item : jsonData.get_child("points_x")) {
        points_x.push_back(item.second.get_value<int>());
    }

    vector<int> points_y;
    for (auto& item : jsonData.get_child("points_y")) {
        points_y.push_back(item.second.get_value<int>());
    }

    vector<int> region_boundary;
    for (auto& item : jsonData.get_child("region_boundary")) {
        region_boundary.push_back(item.second.get_value<int>());
    }

    int num_constraints = jsonData.get<int>("num_constraints");

    vector<vector<int>> additional_constraints;
    for (auto& item : jsonData.get_child("additional_constraints")) {
        vector<int> constraint;
        for (auto& sub_item : item.second) {
            constraint.push_back(sub_item.second.get_value<int>());
        }
        additional_constraints.push_back(constraint);
    }

    vector<Point> points;
    for (int i = 0; i < num_points; ++i) {
        points.push_back(Point(points_x[i], points_y[i]));
    }

    constrained_delaunay_function(points, num_constraints, additional_constraints);

    //make_json();
    return 0;
}