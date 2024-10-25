#include <iostream>
#include <unistd.h>
#include <fcntl.h>
#include <string>
#include <boost/json.hpp>
#include <boost/json/src.hpp>
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
typedef K::Segment_2 Segment;
typedef CDT::Face_handle Face_handle;
typedef CDT::Edge_iterator Edge_iterator;

#define BUFFER 1024

using namespace std;
using namespace boost::json;


int find_point_index(const vector<Point>& points, const Point& target){
    auto it = find(points.begin(), points.end(), target);
    if (it != points.end()){
        return distance(points.begin(), it);
    }
    return -1;
}



int make_json(std::string instance_uid, int num_points, vector<Point> points, vector<Segment> segments){
    std::string content_type = "CG_SHOP_2025_Solution";

    std::stringstream ss;
    ss << "{\n";
    ss << "  \"content_type\": \"" << content_type << "\",\n";
    ss << "  \"instance_uid\": \"" << instance_uid << "\",\n";


    boost::json::array steiner_points_x;
    for (auto it = points.begin() + num_points; it != points.end(); ++it) {

        if (it != points.end()) {
            steiner_points_x.push_back(value((*it)[0]));
        }
    }
    
    ss << "  \"steiner_points_x\": [";
    for (size_t i = 0; i < steiner_points_x.size(); ++i) {
        ss << steiner_points_x[i];
        if (i < steiner_points_x.size() - 1) ss << ", ";
    }
    ss << "],\n";


    boost::json::array steiner_points_y;
    for (auto it = points.begin() + num_points; it != points.end(); ++it){

        if (it != points.end()) {
            steiner_points_y.push_back(value((*it)[1]));
        }
    }
    
    ss << "  \"steiner_points_y\": [";
    for (size_t i = 0; i < steiner_points_y.size(); ++i){
        ss << steiner_points_y[i];
        if (i < steiner_points_y.size() - 1) ss << ", ";
    }
    ss << "],\n";


    ss << "  \"edges\": [\n";
    for (size_t i = 0; i < segments.size(); ++i){
        int index1 = find_point_index(points, segments[i].source());
        int index2 = find_point_index(points, segments[i].target());

        if (index1 != -1 && index2 != -1) {
            ss << "    [" << index1 << ", " << index2 << "]";
            if (i < segments.size() - 1) ss << ",";
            ss << "\n";
        }
    }
    ss << "  ]\n";
    ss << "}\n";

    try{
        std::ofstream output_file("output.json");
        output_file << ss.str();
        output_file.close();
        cout << "output.json created successfully" << endl;
    } catch (const std::exception &e) {
        cerr << "error at creating JSON file: " << e.what() << endl;
        return 1;
    }

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



vector<Segment> get_segments(const CDT& cdt){
    vector<Segment> segments;
    for (Edge_iterator edge = cdt.finite_edges_begin(); edge != cdt.finite_edges_end(); ++edge){
        segments.push_back(cdt.segment(*edge));
    }
    return segments;
}



bool is_obtuse_triangle(const Point& a, const Point& b, const Point& c){
    double ab = CGAL::squared_distance(a, b);
    double bc = CGAL::squared_distance(b, c);
    double ca = CGAL::squared_distance(c, a);

    if ((ab > bc + ca) || (bc > ab + ca) || (ca > ab + bc)) return true;

    return false;
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



Point steiner_at_midpoint(CDT& cdt){
    pair<Point, Point> edge = find_longest_edge(cdt);

    Point p1 = edge.first;
    Point p2 = edge.second;

    Point steiner_point = CGAL::midpoint(p1, p2);
    
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

    char filename[BUFFER] = "";

    cout << "Enter file name:" << endl;
    cin >> filename;

std::ifstream file(filename);
    if (!file.is_open()) {
        cerr << "Error opening file." << endl;
        return 1;
    }
    std::string jsonStr((std::istreambuf_iterator<char>(file)), std::istreambuf_iterator<char>());
    file.close();

    value jsonData;
    try {
        jsonData = parse(jsonStr);
    } catch (const std::exception &e) {
        cerr << "error reading file " << e.what() << endl;
        return 2;
    }

    object& obj = jsonData.as_object();

    std::string instance_uid = obj["instance_uid"].as_string().c_str();
    int num_points = obj["num_points"].as_int64();

    vector<int> points_x;
    for (auto& item : obj["points_x"].as_array()) {
        points_x.push_back(item.as_int64());
    }

    vector<int> points_y;
    for (auto& item : obj["points_y"].as_array()) {
        points_y.push_back(item.as_int64());
    }

    vector<Point> points;
    for (int i = 0; i < num_points; ++i) {
        points.push_back(Point(points_x[i], points_y[i]));
    }

    vector<int> region_boundary;
    for (auto& item : obj["region_boundary"].as_array()) {
        region_boundary.push_back(item.as_int64());
    }

    int num_constraints = obj["num_constraints"].as_int64();

    vector<vector<int>> additional_constraints;
    for (auto& item : obj["additional_constraints"].as_array()) {
        vector<int> constraint;
        for (auto& sub_item : item.as_array()) {
            constraint.push_back(sub_item.as_int64());
        }
        additional_constraints.push_back(constraint);
    }

    CDT cdt = constrained_delaunay_function(points, num_constraints, additional_constraints);

    CGAL::draw(cdt);
    int obtuse_triangle_count = count_obtuse_triangles(cdt);
    cout << "Obtuse triangle count is: " << obtuse_triangle_count << endl;

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
            return 3;
        }
        if (steiner_point[0] != nan("") && steiner_point[1] != nan("")) points.push_back(steiner_point);
        else break;
    }
    CGAL::draw(cdt);
    obtuse_triangle_count = count_obtuse_triangles(cdt);
    cout << "Obtuse triangle count is: " << obtuse_triangle_count << endl;

    vector<Segment> segments = get_segments(cdt);

    make_json(instance_uid, num_points, points, segments);

    return 0;
}