#include <iostream>
#include <unistd.h>
#include <fcntl.h>
#include <string>
#include <boost/json.hpp>
#include <boost/json/src.hpp>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/draw_triangulation_2.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/Polygon_2.h>
#include <cmath>
#include <numeric>

typedef CGAL::Exact_predicates_exact_constructions_kernel EK;
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Triangulation_vertex_base_2<K> Vb;
typedef CGAL::Constrained_triangulation_face_base_2<K> Fb;
typedef CGAL::Triangulation_data_structure_2<Vb,Fb> TDS;
typedef CGAL::Exact_intersections_tag Itag;
typedef CGAL::Constrained_Delaunay_triangulation_2<K, TDS, Itag> CDT;
typedef CGAL::Polygon_2<K> Polygon;
typedef CGAL::Bounded_side Bounded_side;
typedef K::Point_2 Point;
typedef K::Segment_2 Segment;
typedef CDT::Face_handle Face_handle;
typedef CDT::Edge_iterator Edge_iterator;
typedef K::Line_2 Line;
typedef EK::FT FT;

using namespace std;
using namespace boost::json;


//
//TODO PARAMETERS IN JSON
//

std::string gmpz_to_string(const CGAL::Gmpz& value) {
    // Use an ostringstream to convert Gmpz to string
    std::ostringstream oss;
    oss << value;  // CGAL::Gmpz supports stream output
    return oss.str();
}

std::string print_rational(const K::FT& coord) {
    // Convert Lazy_exact_nt to exact (CGAL::Gmpq)
    auto exact_coord = CGAL::exact(coord);  // Evaluates Lazy_exact_nt to CGAL::Gmpq

    // Cast the exact coordinate to CGAL::Gmpq
    const CGAL::Gmpq& gmpq_value = exact_coord;

    // Convert numerator and denominator to strings
    std::string num = gmpz_to_string(gmpq_value.numerator());
    std::string den = gmpz_to_string(gmpq_value.denominator());

    // Construct the string representation
    return num + "/" + den;
}



int find_point_index(const vector<Point>& points, const Point& target){
    //search for the target point in the vector
    auto it = find(points.begin(), points.end(), target);
    //check if the target was found
    if (it != points.end()){
        return distance(points.begin(), it);
    }
    return -1;
}



int make_json(std::string instance_uid, int num_points, vector<Point> points, vector<Segment> segments, int obtuse_triangle_count, std::string method){
    //set the content type in a string for a cleaner code :)
    std::string content_type = "CG_SHOP_2025_Solution";

    //use a stringstream to construct the JSON file contents
    std::stringstream ss;
    ss << "{\n";
    ss << "  \"content_type\": \"" << content_type << "\",\n";
    ss << "  \"instance_uid\": \"" << instance_uid << "\",\n";

    //arrays to store coordinates from Steiner points
    boost::json::array steiner_points_x;
    boost::json::array steiner_points_y;
    for (auto it = points.begin() + num_points; it != points.end(); ++it){ //steiner points are offset by numpoints in points vector
        //push integer Steiner coordinates and check for floating-point
        if ((*it)[0] == (int)(*it)[0] || (*it)[0] == 0.0) steiner_points_x.push_back(value((int)(*it)[0]));
        else steiner_points_x.push_back(value(print_rational((*it)[0])));
        if ((*it)[1] == (int)(*it)[1] || (*it)[1] == 0.0) steiner_points_y.push_back(value((int)(*it)[1]));
        else steiner_points_y.push_back(value(print_rational((*it)[1])));
    }
    
    ss << "  \"steiner_points_x\": [";
    for (size_t i = 0; i < steiner_points_x.size(); ++i){
        ss << steiner_points_x[i];
        if (i < steiner_points_x.size() - 1) ss << ", ";
    }
    ss << "],\n";

    ss << "  \"steiner_points_y\": [";
    for (size_t i = 0; i < steiner_points_y.size(); ++i){
        ss << steiner_points_y[i];
        if (i < steiner_points_y.size() - 1) ss << ", ";
    }
    ss << "],\n";

    ss << "  \"edges\": [\n";
    for (size_t i = 0; i < segments.size(); ++i){
        //find the index of each segment endpoint in the point vector
        int index1 = find_point_index(points, segments[i].source());
        int index2 = find_point_index(points, segments[i].target());

        // check if both points are valid and add them to the stringstream
        if (index1 != -1 && index2 != -1) {
            ss << "    [" << index1 << ", " << index2 << "]";
            if (i < segments.size() - 1) ss << ",";
            ss << "\n";
        }
    }
    ss << "  ]\n";

    ss << "  \"obtuse_count\": \"" << obtuse_triangle_count << "\",\n";
    ss << "  \"method\": \"" << method << "\",\n";

    ss << "}\n";

    //try to write output.json
    try{
        std::ofstream output_file("output.json");
        output_file << ss.str();
        output_file.close();
        cout << "output.json created successfully" << endl;
    } catch (const std::exception &e){
        cerr << "error at creating JSON file: " << e.what() << endl;
        return -1;
    }

    return 0;
}



CDT constrained_delaunay_function(vector<Point> points, const vector<int> region_boundary, vector<vector<int>> additional_constraints, Polygon &polygon){
    //initialize a Constrained Delaunay Triangulation object
    CDT cdt;

    //insert the edges of the region boundary as constraints in the triangulation
    for (size_t i = 0; i < region_boundary.size(); ++i){
        //get the indices of two consecutive points
        int a = region_boundary[i];
        int b = region_boundary[(i + 1) % region_boundary.size()]; //wraps around
        //add the point to the polygon and the segment to the constraints
        polygon.push_back(points[a]);
        cdt.insert_constraint(points[a], points[b]);
    }

    //insert all points in the triangulation
    for (const auto& pt : points){
        cdt.insert(pt);
    }


    //insert the additional constraints in the triangulation
    for (const auto& constraint : additional_constraints){
        int a = constraint[0];
        int b = constraint[1];
        cdt.insert_constraint(points[a], points[b]);
    }

    return cdt;
}



vector<Segment> get_segments(const CDT& cdt){
    //vector to store segments
    vector<Segment> segments;

    //iterate all the edges in the cdt to push them as segments
    for (Edge_iterator edge = cdt.finite_edges_begin(); edge != cdt.finite_edges_end(); ++edge){
        segments.push_back(cdt.segment(*edge));
    }
    return segments;
}



bool is_obtuse_triangle(const Face_handle face, const Polygon polygon){
    //retrieve the three points of the triangle
    Point p1 = face->vertex(0)->point();
    Point p2 = face->vertex(1)->point();
    Point p3 = face->vertex(2)->point();

    //check for points lies outside the boundary
    if (polygon.bounded_side(p1) != CGAL::ON_BOUNDED_SIDE &&
        polygon.bounded_side(p1) != CGAL::ON_BOUNDARY){
        return false;
    }
    if (polygon.bounded_side(p2) != CGAL::ON_BOUNDED_SIDE &&
        polygon.bounded_side(p2) != CGAL::ON_BOUNDARY){
        return false;
    }
    if (polygon.bounded_side(p3) != CGAL::ON_BOUNDED_SIDE &&
        polygon.bounded_side(p3) != CGAL::ON_BOUNDARY){
        return false;
    }

    //find the centroid of the triangle and check if it is within the boundary
    Point centroid = CGAL::centroid(p1, p2, p3);
    if (polygon.bounded_side(centroid) != CGAL::ON_BOUNDED_SIDE){
        return false;
    }

    //calculate the vectors for the triangle's edges
    K::Vector_2 v1 = p2 - p1;
    K::Vector_2 v2 = p3 - p1;
    K::Vector_2 v3 = p3 - p2;

    //calculate the dot products
    double dot1 = v1 * v2;
    double dot2 = -v1 * v3;
    double dot3 = v2 * v3;

    //obtuse angle has negative dot product
    return (dot1 < 0 || dot2 < 0 || dot3 < 0);
}



int count_obtuse_triangles(const CDT& cdt, Polygon polygon){
    int obtuse_triangle_count = 0;
    int triangle_count = 0;

    //iterate all the faces in the cdt
    for (auto face = cdt.finite_faces_begin(); face != cdt.finite_faces_end(); ++face){
        triangle_count++;

        //check if the current triangle is obtuse
        if (is_obtuse_triangle(face, polygon)) obtuse_triangle_count++;
    }

    return obtuse_triangle_count;
}



Face_handle find_obtuse_triangle(CDT& cdt, Polygon polygon, int triangle_count){
    int i = 0;
    //iterate all the faces in the cdt
    for (auto face = cdt.finite_faces_begin(); face != cdt.finite_faces_end(); ++face){

        //check if the current triangle is obtuse
        if (is_obtuse_triangle(face, polygon)){
            if (i == triangle_count) return face;
            i++;
        }
    }
    return nullptr;
}



bool valid_steiner(Point steiner_point, const CDT& cdt){
    //iterate all finite edges in the cdt
    for (auto edge = cdt.finite_edges_begin(); edge != cdt.finite_edges_end(); ++edge){
        //get the source and target points of the edge
        Point source = cdt.segment(*edge).source();
        Point target = cdt.segment(*edge).target();
        if (source == steiner_point || target == steiner_point){
            return false; //if point already exists in the cdt
        }
    }
    return true; //if point does not exist in the cdt
}

Point steiner_at_midpoint(CDT& cdt, Polygon polygon){
    int triangle_count = 0;
    //if triangle results to a non valid point, loops for another triangle
    while(1){
        //find an obtuse triangle in the cdt
        Face_handle triangle = find_obtuse_triangle(cdt, polygon, triangle_count);

        //if no obtuse triangle is found, return null point
        if (triangle == nullptr){
            return Point(nan(""), nan(""));
        }

        //retrieve the three points of the triangle
        Point p1 = triangle->vertex(0)->point();
        Point p2 = triangle->vertex(1)->point();
        Point p3 = triangle->vertex(2)->point();

        //calculate the vectors for the triangle's edges
        K::Vector_2 v1 = p2 - p1;
        K::Vector_2 v2 = p3 - p1;
        K::Vector_2 v3 = p3 - p2;

        //calculate the dot products
        double dot1 = v1 * v2;
        double dot2 = -v1 * v3;
        double dot3 = v2 * v3;

        Point steiner_point;

        //determine which angle is obtuse and calculate the midpoint of the corresponding edge
        if (dot1 < 0) steiner_point = CGAL::midpoint(p2, p3);
        else if (dot2 < 0) steiner_point = CGAL::midpoint(p1, p3);
        else if (dot3 < 0) steiner_point = CGAL::midpoint(p1, p2);

        //check if point already exists there
        if (valid_steiner(steiner_point, cdt)){
            cdt.insert(steiner_point);
            return steiner_point;
        }
        triangle_count++;
    }
}



Point steiner_at_circumcenter(CDT& cdt, Polygon polygon){
    int triangle_count = 0;
    //if triangle results to a non valid point, loops for another triangle
    while(1){
        //find an obtuse triangle in the cdt
        Face_handle triangle = find_obtuse_triangle(cdt, polygon, triangle_count);

        //if no obtuse triangle is found, return null point
        if (triangle == nullptr){
            return Point(nan(""), nan(""));
        }

        //retrieve the three points of the triangle
        Point p1 = triangle->vertex(0)->point();
        Point p2 = triangle->vertex(1)->point();
        Point p3 = triangle->vertex(2)->point();

        Point steiner_point = CGAL::circumcenter(p1, p2, p3);
        //if steiner point is outside of boundary, look for centroid instead
        if (polygon.bounded_side(steiner_point) != CGAL::ON_BOUNDED_SIDE) steiner_point = CGAL::centroid(p1, p2, p3);
        //check if point already exists there
        if (valid_steiner(steiner_point, cdt)){
            cdt.insert(steiner_point);
            return steiner_point;
        }
        triangle_count++;
    }
}



Point steiner_at_projection(CDT& cdt, Polygon polygon){
    int triangle_count = 0;
    //if triangle results to a non valid point, loops for another triangle
    while(1){
        //find an obtuse triangle in the cdt
        Face_handle triangle = find_obtuse_triangle(cdt, polygon, triangle_count);
        //if no obtuse triangle is found, return null point
        if (triangle == nullptr){
            return Point(nan(""), nan(""));
        }

        //retrieve the three points of the triangle
        Point p1 = triangle->vertex(0)->point();
        Point p2 = triangle->vertex(1)->point();
        Point p3 = triangle->vertex(2)->point();

        //calculate the vectors for the triangle's edges
        K::Vector_2 v1 = p2 - p1;
        K::Vector_2 v2 = p3 - p1;
        K::Vector_2 v3 = p3 - p2;

        //calculate the dot products
        double dot1 = v1 * v2;
        double dot2 = -v1 * v3;
        double dot3 = v2 * v3;

        //determine which angle is obtuse and calculate the projection on the corresponding edge
        Point steiner_point;
        if (dot1 < 0){
            Line line(p2, p3);
            steiner_point = line.projection(p1);
        }
        else if (dot2 < 0){
            Line line(p1, p3);
            steiner_point = line.projection(p2);
        }
        else if (dot3 < 0){
            Line line(p1, p2);
            steiner_point = line.projection(p3);
        }

        //check if point already exists there
        if (valid_steiner(steiner_point, cdt)){
            cdt.insert(steiner_point);
            return steiner_point;
        }
        triangle_count++;
    }
}



int main(int argc, char *argv[]){
    //file name insert to open
    std::string filename;
    if (argc > 1) filename = argv[1];
    else {
        cout << "Enter file name:" << endl;
        cin >> filename;
    }

    //open the file and check if it opened successfully
    std::ifstream file(filename);
    if (!file.is_open()){
        cerr << "Error opening file." << endl;
        return 1;
    }

    //read the JSON file as a string
    std::string jsonStr((std::istreambuf_iterator<char>(file)), std::istreambuf_iterator<char>());
    file.close();

    //parse the string into JSON data
    value jsonData;
    try {
        jsonData = parse(jsonStr);
    } catch (const std::exception &e){
        cerr << "error reading file " << e.what() << endl;
        return 2;
    }

    //access the JSON data as object
    object& obj = jsonData.as_object();

    //read the variables
    std::string instance_uid = obj["instance_uid"].as_string().c_str();
    int num_points = obj["num_points"].as_int64();

    vector<int> points_x;
    for (auto& item : obj["points_x"].as_array()){
        points_x.push_back(item.as_int64());
    }

    vector<int> points_y;
    for (auto& item : obj["points_y"].as_array()){
        points_y.push_back(item.as_int64());
    }

    //combine x and y coordinates to create Points
    vector<Point> points;
    for (int i = 0; i < num_points; ++i) {
        points.push_back(Point(points_x[i], points_y[i]));
    }

    vector<int> region_boundary;
    for (auto& item : obj["region_boundary"].as_array()){
        region_boundary.push_back(item.as_int64());
    }

    int num_constraints = obj["num_constraints"].as_int64();

    vector<vector<int>> additional_constraints;
    for (auto& item : obj["additional_constraints"].as_array()){
        vector<int> constraint;
        for (auto& sub_item : item.as_array()) {
            constraint.push_back(sub_item.as_int64());
        }
        additional_constraints.push_back(constraint);
    }

    std::string method = obj["method"].as_string().c_str();
    
    map<std::string, double> parameters;
    for (const auto& param : obj["parameters"].as_object()) {
        //check if the value can be treated as a double
        if (param.value().is_double()) {
            parameters[std::string(param.key())] = param.value().as_double();
        }
        else if (param.value().is_int64()) {
            //convert integers to double
            parameters[std::string(param.key())] = static_cast<double>(param.value().as_int64());
        }
        else {
            //log unsupported types
            cerr << "Warning: Parameter " << param.key() << " is not a numeric type and will be ignored." << endl;
        }
    }
    
    bool delaunay = obj["delaunay"].as_bool();

    //create Constrained Delaunay Triangulation (cdt) and Polygon for region boundary
    Polygon polygon;
    CDT cdt = constrained_delaunay_function(points, region_boundary, additional_constraints, polygon);

    //count obtuse triangles and draw the cdt
    int obtuse_triangle_count = count_obtuse_triangles(cdt, polygon);
    cout << "Obtuse triangle count is: " << obtuse_triangle_count << endl;
    CGAL::draw(cdt);

    //select steiner point placement method
    int smethod;
    if (argc > 2) smethod = atoi(argv[2]);
    else {
        cout << "Enter method number (1,2,3)" << endl;
        cin >> smethod;
    }

    //give maximum amount of allowed steiner points
    int steiner_max;
    if (argc > 3) steiner_max = atoi(argv[3]);
    else {
        cout << "Enter number of max steiner points:" << endl;
        cin >> steiner_max;
    }

    //insert steiner points based on the selected method
    for (int i = 0; i < steiner_max; i++){
        Point steiner_point;
        if (smethod == 1) steiner_point = steiner_at_midpoint(cdt, polygon);
        else if (smethod == 2) steiner_point = steiner_at_circumcenter(cdt, polygon);
        else if (smethod == 3) steiner_point = steiner_at_projection(cdt, polygon);

        //add Steiner point only if it is not null
        if (steiner_point[0] != nan("") && steiner_point[1] != nan("")) points.push_back(steiner_point);
        else break; //break if no valid steiner points are left
    }

    //recount obtuse triangles after inserting steiner points and redraw
    obtuse_triangle_count = count_obtuse_triangles(cdt, polygon);
    cout << "Obtuse triangle count is: " << obtuse_triangle_count << endl;
    cout << points.size() - num_points << " Steiner points added" << endl;
    CGAL::draw(cdt);

    //create segments from cdt for the output file
    vector<Segment> segments = get_segments(cdt);
    make_json(instance_uid, num_points, points, segments, obtuse_triangle_count, method);

    return 0;
}