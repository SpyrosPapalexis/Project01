#include <iostream>
#include <unistd.h>
#include <fcntl.h>
#include <string.h>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>

#define BUFFER 1024

using namespace std;
using namespace boost::property_tree;

int main(void){

    char txtbuffer[BUFFER] = "";
    char filename[BUFFER] = "";

    cout << "Enter file name:" << endl;
    cin >> filename;

    ptree jsonData;
    try {
        read_json(filename, jsonData);
    } catch (const json_parser_error &e) {
        cerr << "error reading file: " << e.what() << endl;
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


    return 0;
}