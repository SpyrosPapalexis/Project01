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
#include <CGAL/Polygon_2.h>
#include <cmath>
#include <numeric>

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
typedef CDT::Vertex_handle Vertex_handle;

using namespace std;
using namespace boost::json;

#define RANDSIZE 1000   //accuracy for random functions

//////////////////////////////////////////////
//TODO                                      //
//change last case in ant colony to polygon //
//enable polygon if we manage to fix it     //
//////////////////////////////////////////////



std::string gmpz_to_string(const CGAL::Gmpz& value){
    //use an ostringstream to convert Gmpz to string
    std::ostringstream oss;
    oss << value;  //CGAL::Gmpz supports stream output
    return oss.str();
}


std::string print_rational(const K::FT& coord){
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


bool is_on_constraint(const Point& source, const Point& target, const Point& steiner_point){
    //check for collinear
    if (CGAL::orientation(source, target, steiner_point) != CGAL::COLLINEAR){
        return false;
    }

    //check if its on the line
    return (std::min(source.x(), target.x()) <= steiner_point.x() &&
            std::max(source.x(), target.x()) >= steiner_point.x() &&
            std::min(source.y(), target.y()) <= steiner_point.y() &&
            std::max(source.y(), target.y()) >= steiner_point.y());
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
        if (cdt.is_constrained(*edge)){
            if (is_on_constraint(source, target, steiner_point)){
                return false; //if point lies on a constrained edge
            }
        }
    }
    return true; //if point does not exist in the cdt
}



Point steiner_at_midpoint(CDT& cdt, Polygon polygon, Face_handle triangle){

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
    return Point(nan(""), nan(""));
}


Point steiner_at_circumcenter(CDT& cdt, Polygon polygon, Face_handle triangle){
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
    return Point(nan(""), nan(""));
}


Point steiner_at_centroid(CDT& cdt, Polygon polygon, Face_handle triangle){
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
    return Point(nan(""), nan(""));
}



bool has_constrained_edges(const CDT& cdt, Face_handle triangle){
    for (int i = 0; i < 3; ++i){
        //create an edge for the face
        CDT::Edge edge(triangle, i);

        //check if edge is constrained
        if (cdt.is_constrained(edge)){
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
    double sum_x = 0;
    double sum_y = 0;
    int num_points = 0;

    // Iterate through the edges of the polygon
    for (const auto& edge : polygon_edges){
        // Get the points of the edge
        Point p1 = edge.first->vertex(edge.second)->point();
        Point p2 = edge.first->vertex((edge.second + 1) % 3)->point();

        // Add the coordinates to the sum
        sum_x += p1.x() + p2.x();
        sum_y += p1.y() + p2.y();

        num_points += 2; // Each edge contributes two points
    }

    // Compute the centroid
    double centroid_x = sum_x / num_points;
    double centroid_y = sum_y / num_points;

    return Point(centroid_x, centroid_y);
}


std::vector<CDT::Edge> find_external_edges(const std::set<Face_handle>& triangles){
    std::set<std::pair<Point, Point>> unique_edges; //normalize the edges
    std::vector<CDT::Edge> result;

    for (const auto& face : triangles){
        for (int i = 0; i < 3; ++i) {
            //get edge points
            Point p1 = face->vertex(i)->point();
            Point p2 = face->vertex((i + 1) % 3)->point();

            //normalize edge
            auto edge = std::make_pair(std::min(p1, p2), std::max(p1, p2));

            //if edge doesnt exist, add it
            if (unique_edges.insert(edge).second){
                result.emplace_back(face, i);
            }
        }
    }

    return result;
}



Point steiner_at_polygon(CDT& cdt, Polygon polygon, Face_handle triangle){
    //dont look for a constrained triangle
    if (has_constrained_edges(cdt, triangle)){
        return Point(nan(""), nan(""));
    }

    cout << "found viable triangle" << endl;

    //find all neighboring obtuse triangles
    std::set<Face_handle> obtuse_triangles = find_neighboring_obtuse_triangles(cdt, triangle, polygon);

    //created needed vectors
    std::vector<CDT::Edge> polygon_edges = find_external_edges(obtuse_triangles);
    cout << "polygon consists of " << polygon_edges.size() / 2 << " edges" << endl;

    //constrain the polygon edges
    std::vector<CDT::Edge> constrained_edges;
    for (const auto& edge : polygon_edges){
        //check if its constrained
        if (!cdt.is_constrained(edge)){
            //if its not constrained, constrain it
            cdt.insert_constraint(edge.first->vertex(edge.second)->point(), edge.first->vertex((edge.second + 1) % 3)->point());
            constrained_edges.push_back(edge); //save newly constrained edge
            cout << "constrain added" << endl;
        }
    }

    cout << "constrained " << constrained_edges.size() << endl;

    //remove points inside the polygon
    //remove points here

    cout << "removed" << endl;

    Point steiner_point = polygon_centroid(polygon_edges);
    
    //check if point already exists there
    bool valid = valid_steiner(steiner_point, cdt);
    if (valid) cdt.insert(steiner_point);

    cout << "steinered" << endl;

    //return removed points
    //return points here

    cout << "returned" << endl;

    //remove constraints
    for (const auto& edge : constrained_edges){
        cdt.remove_constraint(edge.first, edge.second);
        cout << "constrain removed" << endl;
    }

    cout << "unconstrained" << endl;

    if (valid) return steiner_point;
    return Point(nan(""), nan(""));
}



Point steiner_at_projection(CDT& cdt, Polygon polygon, Face_handle triangle){
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
    return Point(nan(""), nan(""));
}


// used to test polygon
// CDT polygon_only(CDT cdt, vector<Point>& points, Polygon polygon, int L){
//     int triangle_count = 0;
//     for(int i = 0; i < L; i++){
//         Face_handle triangle = find_obtuse_triangle(cdt, polygon, triangle_count);
//         if (triangle == nullptr) break;
//         Point steiner_point = steiner_at_polygon(cdt, polygon, triangle);
//         if (!isnan(steiner_point[0]) && !isnan(steiner_point[1])){
//             points.push_back(steiner_point);
//             triangle_count = 0;
//         }else{
//             triangle_count++;
//             i--;
//         }
//     }
//     return cdt;
// }


//brute force method for "delaunay": false
CDT brute_force(CDT cdt, int obtuse_triangle_count, vector<Point>& points, Polygon polygon, int L){
    CDT cdt_best = cdt;
    int best_method = 0;
    int min_obtuse = obtuse_triangle_count;
    vector<Point> method_points[4];
    for (int smethod = 1; smethod <= 4; smethod++){
        CDT cdt_new = cdt;
        method_points[smethod-1] = points;
        int new_obtuse_triangle_count = obtuse_triangle_count;
        int find_obtuse_count = 0;  //counter for find triangle function if steiner point is invalid
        for (int i = 0; i < L && new_obtuse_triangle_count > 0; i++){
            Point steiner_point;
            Face_handle triangle = find_obtuse_triangle(cdt_new, polygon, find_obtuse_count);
            if (triangle == nullptr) break;
            if (smethod == 1) steiner_point = steiner_at_midpoint(cdt_new, polygon, triangle);
            else if (smethod == 2) steiner_point = steiner_at_circumcenter(cdt_new, polygon, triangle);
            else if (smethod == 3) steiner_point = steiner_at_projection(cdt_new, polygon, triangle);
            else if (smethod == 4) steiner_point = steiner_at_centroid(cdt_new, polygon, triangle);
            //else if (smethod == 5) steiner_point = steiner_at_polygon(cdt_new, polygon, triangle);

            //compare the obtuse amount with that of the previous step
            if (count_obtuse_triangles(cdt_new, polygon) > new_obtuse_triangle_count) break;

            new_obtuse_triangle_count = count_obtuse_triangles(cdt_new, polygon);

            //better method found
            if (new_obtuse_triangle_count < min_obtuse){
                best_method = smethod;
                min_obtuse = new_obtuse_triangle_count;
                cdt_best = cdt_new;
            }

            //add Steiner point only if it doesnt overlap
            if (!isnan(steiner_point[0]) && !isnan(steiner_point[1])){
                method_points[smethod-1].push_back(steiner_point);
                find_obtuse_count = 0;
            }else{
                find_obtuse_count++;
                i--;
            }
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
        cdt = cdt_best;
        points = method_points[best_method-1];
    }

    return cdt;
}


//local delaunay search
CDT local_search(CDT cdt, int obtuse_triangle_count, vector<Point>& points, Polygon polygon, int L){
    Point steiner_point;
    Point steiner_point_temp;
    int find_obtuse_count = 0;
    for (int i = 0; i < L && obtuse_triangle_count > 0; i++){
        CDT cdt_best = cdt;
        bool checked = false;

        // check all 5 methods
        CDT cdt_new = cdt;
        Face_handle triangle = find_obtuse_triangle(cdt_new, polygon, find_obtuse_count);
        if (triangle == nullptr) break;
        steiner_point_temp = steiner_at_midpoint(cdt_new, polygon, triangle);
        int new_obtuse_triangle_count = count_obtuse_triangles(cdt_new, polygon);
        if (new_obtuse_triangle_count < obtuse_triangle_count){
            cdt_best = cdt_new;
            steiner_point = steiner_point_temp;
            checked = true;
            obtuse_triangle_count = new_obtuse_triangle_count;
        }

        cdt_new = cdt;
        triangle = find_obtuse_triangle(cdt_new, polygon, find_obtuse_count);
        steiner_point_temp = steiner_at_circumcenter(cdt_new, polygon, triangle);
        new_obtuse_triangle_count = count_obtuse_triangles(cdt_new, polygon);
        if (new_obtuse_triangle_count < obtuse_triangle_count){
            cdt_best = cdt_new;
            steiner_point = steiner_point_temp;
            checked = true;
            obtuse_triangle_count = new_obtuse_triangle_count;
        }

        cdt_new = cdt;
        triangle = find_obtuse_triangle(cdt_new, polygon, find_obtuse_count);
        steiner_point_temp = steiner_at_projection(cdt_new, polygon, triangle);
        new_obtuse_triangle_count = count_obtuse_triangles(cdt_new, polygon);
        if (new_obtuse_triangle_count < obtuse_triangle_count){
            cdt_best = cdt_new;
            steiner_point = steiner_point_temp;
            checked = true;
            obtuse_triangle_count = new_obtuse_triangle_count;
        }

        cdt_new = cdt;
        triangle = find_obtuse_triangle(cdt_new, polygon, find_obtuse_count);
        steiner_point_temp = steiner_at_centroid(cdt_new, polygon, triangle);
        new_obtuse_triangle_count = count_obtuse_triangles(cdt_new, polygon);
        if (new_obtuse_triangle_count < obtuse_triangle_count){
            cdt_best = cdt_new;
            steiner_point = steiner_point_temp;
            checked = true;
            obtuse_triangle_count = new_obtuse_triangle_count;
        }

        //same for 5th unimplimented method

        if (checked == false){
            find_obtuse_count++;
            i--;
            continue;
        }
        cdt = cdt_best;
        obtuse_triangle_count = count_obtuse_triangles(cdt, polygon);
        //add Steiner point only if it doesnt overlap
        if (!isnan(steiner_point[0]) && !isnan(steiner_point[1])){
            points.push_back(steiner_point);
            find_obtuse_count = 0;
        }else{
            find_obtuse_count++;
            i--;
        }
    }
    return cdt;
}


//simulated annealing delaunay method
CDT simulated_annealing(CDT cdt, int obtuse_triangle_count, vector<Point>& points, Polygon polygon, int L, double alpha, double beta){
    Point steiner_point;
    CDT cdt_new = cdt;
    double energy = alpha * obtuse_triangle_count;
    double old_energy;
    double temperature = 1;
    int steiner_count = 0;
    int find_obtuse_count = 0;
    while (temperature >= 0 && obtuse_triangle_count > 0){
        cdt_new = cdt;
        Face_handle triangle = find_obtuse_triangle(cdt_new, polygon, find_obtuse_count);
        if (triangle == nullptr) break;
        old_energy = energy;
        int r = rand()%4;
        if (r == 0) steiner_point = steiner_at_midpoint(cdt_new, polygon, triangle);
        else if (r == 1) steiner_point = steiner_at_circumcenter(cdt_new, polygon, triangle);
        else if (r == 2) steiner_point = steiner_at_projection(cdt_new, polygon, triangle);
        else if (r == 3) steiner_point = steiner_at_centroid(cdt_new, polygon, triangle);
        //else steiner_point = steiner_at_polygon(cdt, polygon, triangle);

        //check for steiner point validity
        if (isnan(steiner_point[0]) && isnan(steiner_point[0])){
            find_obtuse_count++;
            continue;
        }

        obtuse_triangle_count = count_obtuse_triangles(cdt_new, polygon);
        steiner_count++;
        energy = alpha * obtuse_triangle_count + beta * steiner_count;

        if (energy - old_energy < 0){
            cdt = cdt_new;
            points.push_back(steiner_point);
            find_obtuse_count = 0;
        }
        else{
            r = rand()%RANDSIZE;
            if (exp((old_energy - energy)/temperature)*RANDSIZE >= r){
                cdt = cdt_new;
                points.push_back(steiner_point);
                find_obtuse_count = 0;
            }
            else{
                obtuse_triangle_count = count_obtuse_triangles(cdt, polygon);
                //steiner_count--;
                find_obtuse_count++;
            }
        }
        temperature = temperature - 1.0/L;
    }
    return cdt;
}

double radius_to_height(Face_handle triangle){
    //get the vertices of the triangle
    Point p1 = triangle->vertex(0)->point();
    Point p2 = triangle->vertex(1)->point();
    Point p3 = triangle->vertex(2)->point();

    //calculate the lengths of the sides
    double a = sqrt(CGAL::squared_distance(p1, p2));
    double b = sqrt(CGAL::squared_distance(p2, p3));
    double c = sqrt(CGAL::squared_distance(p3, p1));

    //find semi-perimeter to use for Heron's formula
    double s = (a + b + c) / 2.0;

    //find the area using Heron's formula
    double area = sqrt(s * (s - a) * (s - b) * (s - c));

    //find circumradius using area
    double R = (a * b * c) / (4.0 * area);

    //find longest side and height
    double longest_side = max({a, b, c});
    double h = (2.0 * area)/longest_side;

    return R/h;
}


Point improve_triangulation(CDT cdt, Polygon polygon, double xi, double psi, double t[], bool improved[]){
    CDT cdt_new = cdt;
    Point steiner_point;
    //for following arrays:
    // 0 refers to projection
    // 1 refers to circumcenter
    // 2 refers to midpoint
    // 3 refers to mean of adjacent obtuse triangles

    //find a random valid triangle
    Face_handle triangle;
    int obtuse_triangle_count = count_obtuse_triangles(cdt, polygon);
    if (obtuse_triangle_count == 0) return Point(nan(""), nan(""));
    int r;
    do{
        r = rand()%obtuse_triangle_count;
        triangle = find_obtuse_triangle(cdt_new, polygon, r);
    }while (triangle == nullptr);
    double heuristic[4];
    double probability[4];

    //find ratio of current triangle
    double ratio = radius_to_height(triangle);

    //calculate heuristics
    heuristic[0] = max(0.0,(ratio-1)/ratio);
    heuristic[1] = ratio/(2+ratio);
    heuristic[2] = max(0.0,(3-2*ratio)/3);
    heuristic[3] = 1;

    //find sum for the probability
    double sum = 0;
    for (int i = 0; i < 4; i++){
        sum=+ pow(t[i],xi)*pow(heuristic[i],psi);
    }
    
    //find the probability
    for (int sp = 0; sp < 4; sp++){
        probability[sp]=pow(t[sp],xi)*pow(heuristic[sp],psi)/sum*RANDSIZE;
    }

    //find a random method
    r = rand()%RANDSIZE;
    int method;
    if (r <= probability[0]){
        steiner_point = steiner_at_projection(cdt_new, polygon, triangle);
        method = 0;
    }else if (r - probability[0] <= probability[1]){
        steiner_point = steiner_at_circumcenter(cdt_new, polygon, triangle);
        method = 1;
    }else if (r - probability[0] - probability[1] <= probability[2]){
        steiner_point = steiner_at_midpoint(cdt_new, polygon, triangle);
        method = 2;
    }else{
        steiner_point = steiner_at_centroid(cdt_new, polygon, triangle);   //change to polygon later
        method = 3;
    }

    //if the steiner point benefits the cdt, return it. otherwise return nan
    int new_obtuse_triangle_count = count_obtuse_triangles(cdt_new, polygon);
    if (new_obtuse_triangle_count <= obtuse_triangle_count){
        improved[method] = true;
        return steiner_point;
    }
    return Point(nan(""), nan(""));
}


int save_best_triangulation(CDT& cdt, vector<Point>& points, vector<Point> steiner_points, Polygon polygon, double alpha, double beta, int steiner_count){
    CDT cdt_new = cdt;
    int obtuse_triangle_count = count_obtuse_triangles(cdt, polygon);
    double energy = alpha * obtuse_triangle_count + beta * steiner_count;

    for (const Point& p : steiner_points){
        if (valid_steiner(p, cdt_new)){
            cdt_new.insert(p);
        }else{
            continue;
        }
        double old_energy = energy;
        obtuse_triangle_count = count_obtuse_triangles(cdt_new, polygon);
        if (obtuse_triangle_count == 0) energy = 0;
        else energy = alpha * obtuse_triangle_count + beta * steiner_count;

        if (energy < old_energy){
            steiner_count++;
            cdt = cdt_new;
            points.push_back(p);
        }
        cdt_new = cdt;
    }
    return steiner_count;
}


void update_pheromones(int steiner_count, int obtuse_count, bool improved[], double alpha, double beta, double lambda, double t[]){
    for (int i = 0; i < 4; i++){
        double delta;
        if (improved[i]) delta = 1/(1+alpha*obtuse_count+beta*steiner_count);
        else delta = 0;
        t[i] = (1-lambda)*t[i] + delta;
    }
}

//ant colony optimization delaunay method
CDT ant_colony_optimization(CDT cdt, int obtuse_triangle_count, vector<Point>& points, Polygon polygon, int L, double alpha, double beta, int kappa, double lambda, double xi, double psi){
    int obtuse_count = count_obtuse_triangles(cdt, polygon);
    if (obtuse_count == 0) return cdt;
    int steiner_count = 0;
    Point steiner_point;
    double t[4]={1,1,1,1};
    for (int c = 0; c < L; c++){
        
        bool improved[4]={false,false,false,false};
        vector<Point> valid_points;
        
        for (int k = 0; k < kappa; k++){
            Point steiner_point = improve_triangulation(cdt, polygon, xi, psi, t, improved);
            //if non nan point returned from previous function, add it to the valid point vector
            if (!isnan(steiner_point[0]) && !isnan(steiner_point[1])) valid_points.push_back(steiner_point);
        }

        steiner_count = save_best_triangulation(cdt, points, valid_points, polygon, alpha, beta, steiner_count);
        obtuse_count = count_obtuse_triangles(cdt, polygon);
        update_pheromones(steiner_count, obtuse_count, improved, alpha, beta, lambda, t);
    }
    return cdt;
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

    //default parameter values to avoid crashing
    int L = 30;             //number of iterations
    double alpha = 2.0;     //weight of obtuse angles
    double beta = 5.0;      //weight of steiner points
    double xi = 1.0;
    double psi = 3.0;
    double lambda = 0.5;
    int kappa = 10;

    //read parameters from json
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
    
    //delaunay triangulation or brute force
    bool delaunay = false;
    if (obj.contains("delaunay")) delaunay = obj["delaunay"].as_bool();

    //create Constrained Delaunay Triangulation (cdt) and Polygon for region boundary
    Polygon polygon;
    CDT cdt = constrained_delaunay_function(points, region_boundary, additional_constraints, polygon);

    //count obtuse triangles and draw the cdt
    int obtuse_triangle_count = count_obtuse_triangles(cdt, polygon);
    cout << "Obtuse triangle count is: " << obtuse_triangle_count << endl;
    CGAL::draw(cdt);

    if (delaunay == true){
        if (method == "local") cdt = local_search(cdt, obtuse_triangle_count, points, polygon, L); 
        else if (method == "sa") cdt = simulated_annealing(cdt, obtuse_triangle_count, points, polygon, L, alpha, beta);
        else if (method == "ant") cdt = ant_colony_optimization(cdt, obtuse_triangle_count, points, polygon, L, alpha, beta, kappa, lambda, xi, psi);
        else cout << "Invalid method: " << method << endl;
    }
    //brute force method
    else cdt = brute_force(cdt, obtuse_triangle_count, points, polygon, L);

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