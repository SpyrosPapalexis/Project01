#include <iostream>
#include <unistd.h>
#include <fcntl.h>
#include <string.h>

#define BUFFER 1024

using namespace std;

int main(void){

    char buffer[BUFFER] = "";
    char txtbuffer[BUFFER] = "";
    char filename[BUFFER] = "";

    cout << "Enter file name:" << endl;
    cin >> filename;

    strcat(filename,buffer);
    int file;
    if ((file = open(filename, O_RDONLY)) == -1) perror("can't open file");

    if (read(file,txtbuffer,BUFFER) == -1) perror("can't read");

    string instance_uid;
    int numpoints;
    int* points_x;
    int* points_y;
    void* region_boundary;
    int num_constraints;
    void* additional_constraints;





    return 0;
}