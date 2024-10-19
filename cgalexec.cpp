#include <iostream>
#include <unistd.h>
#include <fcntl.h>
#include <string.h>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>

#define BUFFER 1024

using namespace std;
using namespace boost::property_tree;

int make_json(){
    ptree jsonData;

    jsonData.put("content_type", "CG_SHOP_2025_Solution");

    jsonData.put("instance_uid", "unique_instance_id");

    ptree steiner_points_x;
    steiner_points_x.push_back(make_pair("", ptree(std::to_string(123))));  // Ακέραια τιμή ως string
    steiner_points_x.push_back(make_pair("", ptree(std::to_string(789))));  // Ακέραια τιμή ως string
    steiner_points_x.push_back(make_pair("", ptree(std::to_string(1011)))); // Ακέραια τιμή ως string

    jsonData.add_child("steiner_points_x", steiner_points_x);

    // Προσθήκη των steiner_points_y
    ptree steiner_points_y;
    steiner_points_y.push_back(make_pair("", ptree(std::to_string(314))));  // Ακέραια τιμή ως string
    steiner_points_y.push_back(make_pair("", ptree(std::to_string(265))));  // Ακέραια τιμή ως string
    steiner_points_y.push_back(make_pair("", ptree(std::to_string(358))));  // Ακέραια τιμή ως string

    jsonData.add_child("steiner_points_y", steiner_points_y);

    // Προσθήκη των edges (λίστα ακέραιων τιμών)
    ptree edges;
    ptree edge1, edge2, edge3;

    edge1.push_back(make_pair("", ptree(std::to_string(0))));  // Ακέραιες τιμές ως string
    edge1.push_back(make_pair("", ptree(std::to_string(7))));

    edge2.push_back(make_pair("", ptree(std::to_string(7))));
    edge2.push_back(make_pair("", ptree(std::to_string(8))));

    edge3.push_back(make_pair("", ptree(std::to_string(8))));
    edge3.push_back(make_pair("", ptree(std::to_string(9))));

    edges.push_back(make_pair("", edge1));
    edges.push_back(make_pair("", edge2));
    edges.push_back(make_pair("", edge3));

    jsonData.add_child("edges", edges);

    try {
        write_json("output.json", jsonData);
        cout << "output.json created successfully" << endl;
    } catch (const json_parser_error &e) {
        cerr << "couldn't create JSON: " << e.what() << endl;
        return 1;
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

    make_json();
    return 0;
}