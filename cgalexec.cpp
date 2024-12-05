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
typedef CDT::Vertex_handle Vertex_handle;

using namespace std;
using namespace boost::json;

#define RANDSIZE 1000

/////////////////////////////////////
//TODO
//seed for rand in SA
// 
/////////////////////////////////////



std::string gmpz_to_string(const CGAL::Gmpz& value) {
    //use an ostringstream to convert Gmpz to string
    std::ostringstream oss;
    oss << value;  //CGAL::Gmpz supports stream output
    return oss.str();
}


std::string print_rational(const K::FT& coord) {
    //convert Lazy_exact_nt to exact (CGAL::Gmpq)
    auto exact_coord = CGAL::exact(coord);  //evaluates Lazy_exact_nt to CGAL::Gmpq

    //cast the exact coordinate to CGAL::Gmpq
    const CGAL::Gmpq& gmpq_value = exact_coord;

    //convert numerator and denominator to strings
    std::string num = gmpz_to_string(gmpq_value.numerator());
    std::string den = gmpz_to_string(gmpq_value.denominator());

    //construct the string representation
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



int make_json(std::string instance_uid, int num_points, vector<Point> points, vector<Segment> segments, int obtuse_triangle_count, std::string method, map<std::string, double> parameters){
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
    ss << "  ],\n";

    ss << "  \"obtuse_count\": \"" << obtuse_triangle_count << "\",\n";
    ss << "  \"method\": \"" << method << "\",\n";

    //add parameters dynamically
    ss << "  \"parameters\": {\n";
    for (auto it = parameters.begin(); it != parameters.end(); ++it) {
        ss << "    \"" << it->first << "\": " << it->second;
        if (std::next(it) != parameters.end()) ss << ",";
        ss << "\n";
    }
    ss << "  }\n";

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


Point steiner_at_centroid(CDT& cdt, Polygon polygon){
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

        Point steiner_point = CGAL::centroid(p1, p2, p3);
        //check if point already exists there
        if (valid_steiner(steiner_point, cdt)){
            cdt.insert(steiner_point);
            return steiner_point;
        }
        triangle_count++;
    }
}



bool has_constrained_edges(const CDT& cdt, Face_handle triangle){
    for (int i = 0; i < 3; ++i) {
        //create an edge for the face
        CDT::Edge edge(triangle, i);

        //check if edge is constrained
        if (cdt.is_constrained(edge)) {
            return true; //return true if its constrained
        }
    }
    return false; //return false if no edges are constrained
}



std::set<Face_handle> find_neighboring_obtuse_triangles(CDT& cdt, Face_handle triangle, Polygon polygon){
    std::set<Face_handle> obtuse_triangles;
    std::set<Face_handle> visited;
    std::vector<Face_handle> to_visit;

    //start from set triangle
    to_visit.push_back(triangle);

    while (!to_visit.empty()){
        //take last triangle
        Face_handle current = to_visit.back();
        to_visit.pop_back();

        //if triangle is visited, skip it
        if (visited.find(current) != visited.end()){
            continue;
        }

        visited.insert(current);

        //cehck for obtuse
        if (is_obtuse_triangle(current, polygon)){
            obtuse_triangles.insert(current);

            //check for neighbors
            for (int i = 0; i < 3; ++i) {
                Face_handle neighbor = current->neighbor(i);

                //push only if neighbor is obtuse and not visited
                if (!cdt.is_infinite(neighbor) && visited.find(neighbor) == visited.end()){
                    to_visit.push_back(neighbor);
                }
            }
        }
    }

    return obtuse_triangles;
}



Point polygon_centroid(const std::vector<CDT::Edge>& polygon_edges){

    std::vector<Point> polygon_points;

    //get points from edges
    for (const auto& edge : polygon_edges){
        Vertex_handle v1 = edge.first->vertex(edge.second);
        Vertex_handle v2 = edge.first->vertex((edge.second + 1) % 3);
        
        //push back points
        if (std::find(polygon_points.begin(), polygon_points.end(), v1->point()) == polygon_points.end()){
            polygon_points.push_back(v1->point());
        }
        if (std::find(polygon_points.begin(), polygon_points.end(), v2->point()) == polygon_points.end()){
            polygon_points.push_back(v2->point());
        }
    }

    //calculate centroid
    int n = polygon_points.size();
    double A = 0.0;  //area
    double Cx = 0.0; //centroid x
    double Cy = 0.0; //centroid y

    for (int i = 0; i < n; ++i){
        //find each point in a circle
        int j = (i + 1) % n;
        double xi = polygon_points[i].x();
        double yi = polygon_points[i].y();
        double xj = polygon_points[j].x();
        double yj = polygon_points[j].y();

        double cross_product = xi * yj - xj * yi;
        A += cross_product;
        Cx += (xi + xj) * cross_product;
        Cy += (yi + yj) * cross_product;
    }

    A *= 0.5;
    Cx /= (6.0 * A);
    Cy /= (6.0 * A);

    return Point(Cx, Cy);
}



Point steiner_at_polygon(CDT& cdt, Polygon polygon){
    int triangle_count = 0;
    //if triangle results to a non valid point, loops for another triangle
    while(1){
        //find an obtuse triangle in the cdt
        Face_handle triangle = find_obtuse_triangle(cdt, polygon, triangle_count);

        //if no obtuse triangle is found, return null point
        if (triangle == nullptr){
            return Point(nan(""), nan(""));
        }

        if (has_constrained_edges(cdt, triangle)){
            triangle_count++;
            continue;
        }

        cout << "found viable triangle" << endl;

        //find all neighboring obtuse triangles
        std::set<Face_handle> obtuse_triangles = find_neighboring_obtuse_triangles(cdt, triangle, polygon);

        //created needed vectors
        std::vector<CDT::Edge> polygon_edges;

        int k = 0;
        for (const auto& tri : obtuse_triangles){
            for (int i = 0; i < 3; ++i){
                //extract edges from triangles
                CDT::Edge edge(tri, i);
                polygon_edges.push_back(edge);
            }
            k++;
        }
        cout << "polygon consists of " << k << " triangles" << endl;

        //constrain the polygon edges
        std::vector<CDT::Edge> constrained_edges;
        for (const auto& edge : polygon_edges){
            //constrain each edge
            cdt.insert_constraint(edge.first->vertex(edge.second)->point(), edge.first->vertex((edge.second + 1) % 3)->point());
            constrained_edges.push_back(edge); //save constrained edges
        }

        //remove inner edges while saving them
        std::vector<CDT::Edge> removed_edges;
        for (auto& edge : cdt.finite_edges()){
            //check if edge is constrained
            if (!cdt.is_constrained(edge)){
                removed_edges.push_back(edge);
                cdt.remove_constraint(edge.first, edge.second);
            }
        }

        Point steiner_point = polygon_centroid(polygon_edges);
        
        //check if point already exists there
        if (valid_steiner(steiner_point, cdt)){

            cdt.insert(steiner_point);

            //return removed edges
            for (const auto& edge : removed_edges) cdt.insert_constraint(edge.first->vertex(edge.second)->point(), edge.first->vertex((edge.second + 1) % 3)->point());

            //remove constraints
            for (const auto& edge : constrained_edges) cdt.remove_constraint(edge.first, edge.second);

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
    srand(time(NULL));

    //file name insert to open
    std::string filename;
    if (argc > 1) filename = argv[1];
    else{
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
    try{
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
        for (auto& sub_item : item.as_array()){
            constraint.push_back(sub_item.as_int64());
        }
        additional_constraints.push_back(constraint);
    }

    //desired delaunay method
    std::string method = "";
    if (obj.contains("method")) method = obj["method"].as_string().c_str();

    //default parameters to avoid crashing
    int L = 1;              //number of iterations
    double alpha = 0.0;     //weight of obtuse angles
    double beta = 0.0;      //weight of steiner points
    double kappa = 0.0;
    double lambda = 0.0;
    double psi = 0.0;
    double xi = 0.0;

    //desired parameters
    map<std::string, double> parameters;
    if (obj.contains("parameters")) {
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
            if (param.key() == "L") L = parameters[std::string(param.key())];
            if (param.key() == "alpha") alpha = parameters[std::string(param.key())];
            if (param.key() == "beta") beta = parameters[std::string(param.key())];
            if (param.key() == "kappa") kappa = parameters[std::string(param.key())];
            if (param.key() == "lambda") lambda = parameters[std::string(param.key())];
            if (param.key() == "psi") psi = parameters[std::string(param.key())];
            if (param.key() == "xi") xi = parameters[std::string(param.key())];
        }
    }
    
    //use delaunay method or not
    bool delaunay = false;
    if (obj.contains("delaunay")) delaunay = obj["delaunay"].as_bool();

    //create Constrained Delaunay Triangulation (cdt) and Polygon for region boundary
    Polygon polygon;
    CDT cdt = constrained_delaunay_function(points, region_boundary, additional_constraints, polygon);

    //count obtuse triangles and draw the cdt
    int obtuse_triangle_count = count_obtuse_triangles(cdt, polygon);
    cout << "Obtuse triangle count is: " << obtuse_triangle_count << endl;
    CGAL::draw(cdt);

    //give maximum amount of allowed steiner points

    int new_obtuse_triangle_count;
    if (delaunay == true) {
        //local search
        if (method == "local"){
            for (int i = 0; i < L && obtuse_triangle_count > 0; i++){

                Point steiner_point;
                Point steiner_point_temp;

                CDT cdt_best = cdt;
                
                bool checked = false;

                // check all 5 methods
                CDT cdt_new = cdt;
                steiner_point_temp = steiner_at_midpoint(cdt_new, polygon);
                new_obtuse_triangle_count = count_obtuse_triangles(cdt_new, polygon);
                if (new_obtuse_triangle_count <= obtuse_triangle_count){
                    cdt_best = cdt_new;
                    steiner_point = steiner_point_temp;
                    checked = true;
                }

                cdt_new = cdt;
                steiner_point_temp = steiner_at_circumcenter(cdt_new, polygon);
                new_obtuse_triangle_count = count_obtuse_triangles(cdt_new, polygon);
                if (new_obtuse_triangle_count <= obtuse_triangle_count){
                    cdt_best = cdt_new;
                    steiner_point = steiner_point_temp;
                    checked = true;
                }

                cdt_new = cdt;
                steiner_point_temp = steiner_at_projection(cdt_new, polygon);
                new_obtuse_triangle_count = count_obtuse_triangles(cdt_new, polygon);
                if (new_obtuse_triangle_count <= obtuse_triangle_count){
                    cdt_best = cdt_new;
                    steiner_point = steiner_point_temp;
                    checked = true;
                }

                cdt_new = cdt;
                steiner_point_temp = steiner_at_centroid(cdt_new, polygon);
                new_obtuse_triangle_count = count_obtuse_triangles(cdt_new, polygon);
                if (new_obtuse_triangle_count <= obtuse_triangle_count){
                    cdt_best = cdt_new;
                    steiner_point = steiner_point_temp;
                    checked = true;
                }

                //same for 5th unimplimented method

                if (checked == false) break;
                cdt = cdt_best;
                obtuse_triangle_count = count_obtuse_triangles(cdt, polygon);
                //add Steiner point only if it is not null
                if (steiner_point[0] != nan("") && steiner_point[1] != nan("")) points.push_back(steiner_point);
                else break; //break if no valid steiner points are left
            }
        }
        //simulated annealing
        else if (method == "sa") {
            Point steiner_point;
            CDT cdt_new = cdt;
            double epsilon = alpha * obtuse_triangle_count;
            double old_epsilon;
            double temperature = 1;
            int steiner_count = 0;
            while (temperature >= 0 && obtuse_triangle_count > 0){
                old_epsilon = epsilon;
                int r = rand()%4;
                if (r == 0){
                    steiner_point = steiner_at_midpoint(cdt_new, polygon);
                }else if (r == 1){
                    steiner_point = steiner_at_circumcenter(cdt_new, polygon);
                }else if (r == 2){
                    steiner_point = steiner_at_projection(cdt_new, polygon);
                }else if (r == 3){
                    steiner_point = steiner_at_centroid(cdt_new, polygon);
                // }else{
                //     steiner_point = steiner_at_(cdt, polygon);
                }
                obtuse_triangle_count = count_obtuse_triangles(cdt_new, polygon);
                steiner_count++;
                epsilon = alpha * obtuse_triangle_count + beta * steiner_count;

                if (epsilon - old_epsilon < 0){
                    cdt = cdt_new;
                    points.push_back(steiner_point);
                }else{
                    r = rand()%RANDSIZE;
                    if (exp((old_epsilon - epsilon)/temperature)*RANDSIZE >= r){
                        cdt = cdt_new;
                        points.push_back(steiner_point);
                    }
                    else{
                        obtuse_triangle_count = count_obtuse_triangles(cdt, polygon);
                        cdt_new = cdt;
                    }
                }
                temperature = temperature - 1.0/L;
            }
        }
        //ant colony
        else if (method == "ant") {
            cout << method << endl;
        }
        else cout << "Invalid method: " << method << endl;
    }
    //non-delaunay methods
    else {
        CDT cdt_best = cdt;
        int best_method = 0;
        int min_obtuse = obtuse_triangle_count;
        vector<Point> method_points[5];
        for (int smethod = 1; smethod <= 5; smethod++){
            CDT cdt_new = cdt;
            method_points[smethod-1] = points;
            new_obtuse_triangle_count = obtuse_triangle_count;
            for (int i = 0; i < L && new_obtuse_triangle_count > 0; i++){
                Point steiner_point;
                if (smethod == 1) steiner_point = steiner_at_midpoint(cdt_new, polygon);
                else if (smethod == 2) steiner_point = steiner_at_circumcenter(cdt_new, polygon);
                else if (smethod == 3) steiner_point = steiner_at_projection(cdt_new, polygon);
                else if (smethod == 4) steiner_point = steiner_at_centroid(cdt_new, polygon);
                else if (smethod == 5) steiner_point = steiner_at_polygon(cdt_new, polygon);

                //compare the obtuse amount with that of the previous step
                if (count_obtuse_triangles(cdt_new, polygon) > new_obtuse_triangle_count) break;

                new_obtuse_triangle_count = count_obtuse_triangles(cdt_new, polygon);

                //better method found
                if (new_obtuse_triangle_count < min_obtuse){
                    best_method = smethod;
                    min_obtuse = new_obtuse_triangle_count;
                    cdt_best = cdt_new;
                }

                //add Steiner point only if it is not null
                if (steiner_point[0] != nan("") && steiner_point[1] != nan("")) method_points[smethod-1].push_back(steiner_point);
                else break; //break if no valid steiner points are left
            }
        }

        cout << "Best method is: ";
        if (best_method == 0) cout << "None";
        else if (best_method == 1) cout << "Steiner at midpoint";
        else if (best_method == 2) cout << "Steiner at circumcenter";
        else if (best_method == 3) cout << "Steiner at projection";
        else if (best_method == 4) cout << "Steiner at centroid";
        cout << endl;
        
        //return best method
        if (best_method != 0){
            obtuse_triangle_count = min_obtuse;
            cdt = cdt_best;
            points = method_points[best_method-1];
        }

        // for (int i = 0; i < L && obtuse_triangle_count > 0; i++){
        //     steiner_point = steiner_at_polygon(cdt, polygon);

        //     //add Steiner point only if it is not null
        //     if (steiner_point[0] != nan("") && steiner_point[1] != nan("")) points.push_back(steiner_point);
        //     else break; //break if no valid steiner points are left
        // }
    }

    //recount obtuse triangles after inserting steiner points and redraw
    obtuse_triangle_count = count_obtuse_triangles(cdt, polygon);
    cout << "Obtuse triangle count is: " << obtuse_triangle_count << endl;
    cout << "Steiner points added: " << points.size() - num_points << endl;
    CGAL::draw(cdt);

    //create segments from cdt for the output file
    vector<Segment> segments = get_segments(cdt);
    make_json(instance_uid, num_points, points, segments, obtuse_triangle_count, method, parameters);

    return 0;
}