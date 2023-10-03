#include<mpi.h>
#include<vector>
#include<iostream>
#include"shapefunction.hpp"
#include<fstream>
#include<time.h>

using namespace std;

void set_field(vector<vector<double>> &node, vector<vector<int>> &element, int n, double dx)
{
    node.resize(n*n);
    for(int i=0; i<node.size(); i++){
        node[i].resize(2);
    }
    for(int i=0; i<n; i++){
        for(int j=0; j<n; j++){
            node[i*n+j][0]= j*dx;
            node[i*n+j][1]= i*dx;
        }
    }
    element.resize((n-1)*(n-1));
    for(int i=0; i<n-1; i++){
        for(int j=0; j<n-1; j++){
            element[i*(n-1)+j].push_back(i*n+j);
            element[i*(n-1)+j].push_back(i*n+j+1);
            element[i*(n-1)+j].push_back(i*n+j+(n+1));
            element[i*(n-1)+j].push_back(i*n+j+n);   
        }
    }
}

void set_global_node_number(vector<int> &local_data, int local_grid_x, int local_grid_y, int local_size_x, int local_size_y, int corner, int point_num, int Gx, int Gy)
{
    //case0
    if(local_grid_x==0 && local_grid_y==0){
        for(int i=0; i<local_size_y+1; i++){
            for(int j=0; j<local_size_x+1; j++){
                if(i==0 && j==0){
                    local_data.push_back(corner);
                }
                else if(i!=local_size_y && j==0){
                    local_data.push_back(corner+point_num*i);
                }
                else if(j==local_size_x){
                    local_data.push_back(0);
                }
                else if(i==local_size_y){
                    local_data.push_back(0);
                }
                else{
                    local_data.push_back(local_data[local_data.size()-1]+1);
                }
            }
        }
    }
    //case2
    else if(local_grid_x==Gx-1 && local_grid_y==0){
        for(int i=0; i<local_size_y+1; i++){
            for(int j=0; j<local_size_x+1; j++){
                if(j==0){
                    local_data.push_back(0);
                }
                else if(i==local_size_y){
                    local_data.push_back(0);
                }
                else if(i==0 && j==1){
                    local_data.push_back(corner);
                }
                else if(i!=0 && j==1){
                    local_data.push_back(corner+point_num*i);
                }
                else{
                    local_data.push_back(local_data[local_data.size()-1]+1);
                }
            }
        }
    }
    //case1
    else if(local_grid_y==0 && local_grid_x<Gx-1 && local_grid_x>0){
        for(int i=0; i<local_size_y+1; i++){
            for(int j=0; j<local_size_x+2; j++){
                if(j==0){
                    local_data.push_back(0);
                }
                else if(i==local_size_y){
                    local_data.push_back(0);
                }
                else if(j==local_size_x+1){
                    local_data.push_back(0);
                }
                else if(i==0 && j==1){
                    local_data.push_back(corner);
                }
                else if(i!=0 && j==1){
                    local_data.push_back(corner+point_num*i);
                }
                else{
                    local_data.push_back(local_data[local_data.size()-1]+1);
                }
            }
        }
    }
    //case6
    else if(local_grid_x==0 && local_grid_y==Gy-1){
        for(int i=0; i<local_size_y+1; i++){
            for(int j=0; j<local_size_x+1; j++){
                if(i==0){
                    local_data.push_back(0);
                }
                else if(j==local_size_x){
                    local_data.push_back(0);
                }
                else if(i==1 && j==0){
                    local_data.push_back(corner);
                }
                else if(i!=1 && j==0){
                    local_data.push_back(corner+point_num*(i-1));
                }
                else{
                    local_data.push_back(local_data[local_data.size()-1]+1);
                }
            }
        }
    }
    //case 3
    else if(local_grid_x==0 && local_grid_y<Gy-1 && local_grid_y>0){
        for(int i=0; i<local_size_y+2; i++){
            for(int j=0; j<local_size_x+1; j++){
                if(i==0){
                    local_data.push_back(0);
                }
                else if(i==local_size_y+1){
                    local_data.push_back(0);
                }
                else if(j==local_size_x){
                    local_data.push_back(0);
                }
                else if(i==1 && j==0){
                    local_data.push_back(corner);
                }
                else if(i!=1 && j==0){
                    local_data.push_back(corner+point_num*(i-1));
                }
                else{
                    local_data.push_back(local_data[local_data.size()-1]+1);
                }
            }
        }
    }
    //case 8
    else if(local_grid_x==Gx-1 && local_grid_y==Gy-1){
        for(int i=0; i<local_size_y+1; i++){
            for(int j=0; j<local_size_x+1; j++){
                if(i==0){
                    local_data.push_back(0);
                }
                else if(j==0){
                    local_data.push_back(0);
                }
                else if(i==1 && j==1){
                    local_data.push_back(corner);
                }
                else if(i!=1 && j==1){
                    local_data.push_back(corner+point_num*(i-1));
                }
                else{
                    local_data.push_back(local_data[local_data.size()-1]+1);
                }
            }
        }
    }
    //case7
    else if(local_grid_x<Gx-1 && local_grid_x>0 && local_grid_y==Gy-1){
        for(int i=0; i<local_size_y+1; i++){
            for(int j=0; j<local_size_x+2; j++){
                if(i==0){
                    local_data.push_back(0);
                }
                else if(j==0){
                    local_data.push_back(0);
                }
                else if(j==local_size_x+1){
                    local_data.push_back(0);
                }
                else if(i==1 && j==1){
                    local_data.push_back(corner);
                }
                else if(i!=1 && j==1){
                    local_data.push_back(corner+point_num*(i-1));
                }
                else{
                    local_data.push_back(local_data[local_data.size()-1]+1);
                }
            }
        }
    }
    //case5
    else if(local_grid_x==Gx-1 && local_grid_y>0 && local_grid_y<Gy-1){
        for(int i=0; i<local_size_y+1; i++){
            for(int j=0; j<local_size_x+2; j++){
                if(i==0){
                    local_data.push_back(0);
                }
                else if(j==0){
                    local_data.push_back(0);
                }
                else if(i==local_size_y+1){
                    local_data.push_back(0);
                }
                else if(i==1 && j==1){
                    local_data.push_back(corner);
                }
                else if(i!=1 && j==1){
                    local_data.push_back(corner+point_num*(i-1));
                }
                else{
                    local_data.push_back(local_data[local_data.size()-1]+1);
                }
            }
        }
    }
    //case4
    else if(local_grid_x<Gx-1 && local_grid_x>0 && local_grid_y>0 && local_grid_y<Gy-1){
        for(int i=0; i<local_size_y+2; i++){
            for(int j=0; j<local_size_x+2; j++){
                if(i==0){
                    local_data.push_back(0);
                }
                else if(j==0){
                    local_data.push_back(0);
                }
                else if(j==local_size_x+1){
                    local_data.push_back(0);
                }
                else if(i==local_size_y+1){
                    local_data.push_back(0);
                }
                else if(i==1 && j==1){
                    local_data.push_back(corner);
                }
                else if(i!=1 && j==1){
                    local_data.push_back(corner+point_num*(i-1));
                }
                else{
                    local_data.push_back(local_data[local_data.size()-1]+1);
                }
            }
        }
    }
}

void send_recv_data_x_int(vector<int> &local_data, int rank, int local_grid_x, int local_grid_y, int local_size_x, int local_size_y, int Gx, int Gy)
{
    //のりしろの通信
    MPI_Status st;
    //case0
    if(local_grid_x==0 && local_grid_y==0){
        vector<int> send_test_right;
        vector<int> recv_test_right(local_size_y);
        for (int i = 0; i < local_size_y; i++) {
            send_test_right.push_back(local_data[(local_size_x+1)*i+(local_size_x-1)]);
        }
        MPI_Sendrecv(send_test_right.data(), local_size_y, MPI_INT, rank+1, 0, recv_test_right.data(), local_size_y, MPI_INT, rank+1, 0, MPI_COMM_WORLD, &st);
        for(int i=0; i<recv_test_right.size(); i++){
            local_data[(local_size_x+1)*i+local_size_x]=recv_test_right[i];
        }
    }
    //case1
    if(local_grid_y==0 && local_grid_x<Gx-1 && local_grid_x>0){
        vector<int> send_test_left;
        vector<int> send_test_right;
        vector<int> recv_test_left(local_size_y);
        vector<int> recv_test_right(local_size_y);
        for (int i = 0; i < local_size_y; i++) {
            send_test_left.push_back(local_data[(local_size_x+2)*i+1]);
            send_test_right.push_back(local_data[(local_size_x+2)*i+local_size_x]);
        }
        MPI_Sendrecv(send_test_left.data(), local_size_y, MPI_INT, rank-1, 0, recv_test_left.data(), local_size_y, MPI_INT, rank-1, 0, MPI_COMM_WORLD, &st);
        MPI_Sendrecv(send_test_right.data(), local_size_y, MPI_INT, rank+1, 0, recv_test_right.data(), local_size_y, MPI_INT, rank+1, 0, MPI_COMM_WORLD, &st);
        for(int i=0; i<recv_test_right.size(); i++){
            local_data[(local_size_x+2)*i+local_size_x+1]=recv_test_right[i];
        }
        for(int i=0; i<recv_test_left.size(); i++){
            local_data[(local_size_x+2)*i]=recv_test_left[i];
        }
    }
    //case2
    if(local_grid_x==Gx-1 && local_grid_y==0){
        vector<int> send_test_left;
        vector<int> recv_test_left(local_size_y);
        for (int i = 0; i < local_size_y; i++) {
            send_test_left.push_back(local_data[(local_size_x+1)*i+1]);
        }
        MPI_Sendrecv(send_test_left.data(), local_size_y, MPI_INT, rank-1, 0, recv_test_left.data(), local_size_y, MPI_INT, rank-1, 0, MPI_COMM_WORLD, &st);
        for(int i=0; i<recv_test_left.size(); i++){
            local_data[(local_size_x+1)*i]=recv_test_left[i];
        }
    }
    //case6
    if(local_grid_x==0 && local_grid_y==Gy-1){
        vector<int> send_test_right;
        vector<int> recv_test_right(local_size_y);
        for (int i = 0; i < local_size_y; i++) {
            send_test_right.push_back(local_data[(local_size_x+1)*(i+1)+(local_size_x-1)]);
        }
        MPI_Sendrecv(send_test_right.data(), local_size_y, MPI_INT, rank+1, 0, recv_test_right.data(), local_size_y, MPI_INT, rank+1, 0, MPI_COMM_WORLD, &st);
        for(int i=0; i<recv_test_right.size(); i++){
            local_data[(local_size_x+1)*(i+1)+local_size_x]=recv_test_right[i];
        }
    }
    //case7
    if(local_grid_x<Gx-1 && local_grid_x>0 && local_grid_y==Gy-1){
        vector<int> send_test_left;
        vector<int> send_test_right;
        vector<int> recv_test_left(local_size_y);
        vector<int> recv_test_right(local_size_y);
        for (int i = 0; i < local_size_y; i++) {
            send_test_left.push_back(local_data[(local_size_x+2)*(i+1)+1]);
            send_test_right.push_back(local_data[(local_size_x+2)*(i+1)+local_size_x]);
        }
        MPI_Sendrecv(send_test_left.data(), local_size_y, MPI_INT, rank-1, 0, recv_test_left.data(), local_size_y, MPI_INT, rank-1, 0, MPI_COMM_WORLD, &st);
        MPI_Sendrecv(send_test_right.data(), local_size_y, MPI_INT, rank+1, 0, recv_test_right.data(), local_size_y, MPI_INT, rank+1, 0, MPI_COMM_WORLD, &st);
        for(int i=0; i<recv_test_right.size(); i++){
            local_data[(local_size_x+2)*(i+1)+local_size_x+1]=recv_test_right[i];
        }
        for(int i=0; i<recv_test_left.size(); i++){
            local_data[(local_size_x+2)*(i+1)]=recv_test_left[i];
        }
    }
    //case8
    if(local_grid_x==Gx-1 && local_grid_y==Gy-1){
        vector<int> send_test_left;
        vector<int> recv_test_left(local_size_y);
        for (int i = 0; i < local_size_y; i++) {
            send_test_left.push_back(local_data[(local_size_x+1)*(i+1)+1]);
        }
        MPI_Sendrecv(send_test_left.data(), local_size_y, MPI_INT, rank-1, 0, recv_test_left.data(), local_size_y, MPI_INT, rank-1, 0, MPI_COMM_WORLD, &st);
        for(int i=0; i<recv_test_left.size(); i++){
            local_data[(local_size_x+1)*(i+1)]=recv_test_left[i];
        }
    }    
}

void send_recv_data_x_double(vector<double> &local_data, int rank, int local_grid_x, int local_grid_y, int local_size_x, int local_size_y, int Gx, int Gy)
{
    //のりしろの通信
    MPI_Status st;
    //case0
    if(local_grid_x==0 && local_grid_y==0){
        vector<double> send_test_right;
        vector<double> recv_test_right(local_size_y);
        for (int i = 0; i < local_size_y; i++) {
            send_test_right.push_back(local_data[(local_size_x+1)*i+(local_size_x-1)]);
        }
        MPI_Sendrecv(send_test_right.data(), local_size_y, MPI_DOUBLE, rank+1, 0, recv_test_right.data(), local_size_y, MPI_DOUBLE, rank+1, 0, MPI_COMM_WORLD, &st);
        for(int i=0; i<recv_test_right.size(); i++){
            local_data[(local_size_x+1)*i+local_size_x]=recv_test_right[i];
        }
    }
    //case1
    if(local_grid_y==0 && local_grid_x<Gx-1 && local_grid_x>0){
        vector<double> send_test_left;
        vector<double> send_test_right;
        vector<double> recv_test_left(local_size_y);
        vector<double> recv_test_right(local_size_y);
        for (int i = 0; i < local_size_y; i++) {
            send_test_left.push_back(local_data[(local_size_x+2)*i+1]);
            send_test_right.push_back(local_data[(local_size_x+2)*i+local_size_x]);
        }
        MPI_Sendrecv(send_test_left.data(), local_size_y, MPI_DOUBLE, rank-1, 0, recv_test_left.data(), local_size_y, MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD, &st);
        MPI_Sendrecv(send_test_right.data(), local_size_y, MPI_DOUBLE, rank+1, 0, recv_test_right.data(), local_size_y, MPI_DOUBLE, rank+1, 0, MPI_COMM_WORLD, &st);
        for(int i=0; i<recv_test_right.size(); i++){
            local_data[(local_size_x+2)*i+local_size_x+1]=recv_test_right[i];
        }
        for(int i=0; i<recv_test_left.size(); i++){
            local_data[(local_size_x+2)*i]=recv_test_left[i];
        }
    }
    //case2
    if(local_grid_x==Gx-1 && local_grid_y==0){
        vector<double> send_test_left;
        vector<double> recv_test_left(local_size_y);
        for (int i = 0; i < local_size_y; i++) {
            send_test_left.push_back(local_data[(local_size_x+1)*i+1]);
        }
        MPI_Sendrecv(send_test_left.data(), local_size_y, MPI_DOUBLE, rank-1, 0, recv_test_left.data(), local_size_y, MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD, &st);
        for(int i=0; i<recv_test_left.size(); i++){
            local_data[(local_size_x+1)*i]=recv_test_left[i];
        }
    }
    //case6
    if(local_grid_x==0 && local_grid_y==Gy-1){
        vector<double> send_test_right;
        vector<double> recv_test_right(local_size_y);
        for (int i = 0; i < local_size_y; i++) {
            send_test_right.push_back(local_data[(local_size_x+1)*(i+1)+(local_size_x-1)]);
        }
        MPI_Sendrecv(send_test_right.data(), local_size_y, MPI_DOUBLE, rank+1, 0, recv_test_right.data(), local_size_y, MPI_DOUBLE, rank+1, 0, MPI_COMM_WORLD, &st);
        for(int i=0; i<recv_test_right.size(); i++){
            local_data[(local_size_x+1)*(i+1)+local_size_x]=recv_test_right[i];
        }
    }
    //case7
    if(local_grid_x<Gx-1 && local_grid_x>0 && local_grid_y==Gy-1){
        vector<double> send_test_left;
        vector<double> send_test_right;
        vector<double> recv_test_left(local_size_y);
        vector<double> recv_test_right(local_size_y);
        for (int i = 0; i < local_size_y; i++) {
            send_test_left.push_back(local_data[(local_size_x+2)*(i+1)+1]);
            send_test_right.push_back(local_data[(local_size_x+2)*(i+1)+local_size_x]);
        }
        MPI_Sendrecv(send_test_left.data(), local_size_y, MPI_DOUBLE, rank-1, 0, recv_test_left.data(), local_size_y, MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD, &st);
        MPI_Sendrecv(send_test_right.data(), local_size_y, MPI_DOUBLE, rank+1, 0, recv_test_right.data(), local_size_y, MPI_DOUBLE, rank+1, 0, MPI_COMM_WORLD, &st);
        for(int i=0; i<recv_test_right.size(); i++){
            local_data[(local_size_x+2)*(i+1)+local_size_x+1]=recv_test_right[i];
        }
        for(int i=0; i<recv_test_left.size(); i++){
            local_data[(local_size_x+2)*(i+1)]=recv_test_left[i];
        }
    }
    //case8
    if(local_grid_x==Gx-1 && local_grid_y==Gy-1){
        vector<double> send_test_left;
        vector<double> recv_test_left(local_size_y);
        for (int i = 0; i < local_size_y; i++) {
            send_test_left.push_back(local_data[(local_size_x+1)*(i+1)+1]);
        }
        MPI_Sendrecv(send_test_left.data(), local_size_y, MPI_DOUBLE, rank-1, 0, recv_test_left.data(), local_size_y, MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD, &st);
        for(int i=0; i<recv_test_left.size(); i++){
            local_data[(local_size_x+1)*(i+1)]=recv_test_left[i];
        }
    }    
}

void send_recv_data_y_int(vector<int> &local_data, int rank, int local_grid_x, int local_grid_y, int local_size_x, int local_size_y, int Gx, int Gy)
{
    //のりしろの通信
    MPI_Status st;
    //case0
    if(local_grid_x==0 && local_grid_y==0){
        vector<int> send_test_lower;
        vector<int> recv_test_lower(local_size_x+1);
        for (int i = 0; i < local_size_x+1; i++) {
            send_test_lower.push_back(local_data[(local_size_x+1)*(local_size_y-1)+i]);
        }
        MPI_Sendrecv(send_test_lower.data(), local_size_x+1, MPI_INT, rank+Gx, 0, recv_test_lower.data(), local_size_x+1, MPI_INT, rank+Gx, 0, MPI_COMM_WORLD, &st);
        for(int i=0; i<recv_test_lower.size(); i++){
            local_data[(local_size_x+1)*local_size_y+i]=recv_test_lower[i];
        }
    }
    //case1
    if(local_grid_y==0 && local_grid_x<Gx-1 && local_grid_x>0){
        vector<int> send_test_lower;
        vector<int> recv_test_lower(local_size_x+2);
        for (int i = 0; i < local_size_x+2; i++) {
            send_test_lower.push_back(local_data[(local_size_x+2)*(local_size_y-1)+i]);
        }
        MPI_Sendrecv(send_test_lower.data(), local_size_x+2, MPI_INT, rank+Gx, 0, recv_test_lower.data(), local_size_x+2, MPI_INT, rank+Gx, 0, MPI_COMM_WORLD, &st);
        for(int i=0; i<recv_test_lower.size(); i++){
            local_data[(local_size_x+2)*local_size_y+i]=recv_test_lower[i];
        }
    }
    //case2
    if(local_grid_x==Gx-1 && local_grid_y==0){
        vector<int> send_test_lower;
        vector<int> recv_test_lower(local_size_x+1);
        for (int i = 0; i < local_size_x+1; i++) {
            send_test_lower.push_back(local_data[(local_size_x+1)*(local_size_y-1)+i]);
        }
        MPI_Sendrecv(send_test_lower.data(), local_size_x+1, MPI_INT, rank+Gx, 0, recv_test_lower.data(), local_size_x+1, MPI_INT, rank+Gx, 0, MPI_COMM_WORLD, &st);
        for(int i=0; i<recv_test_lower.size(); i++){
            local_data[(local_size_x+1)*local_size_y+i]=recv_test_lower[i];
        }
    }
    //case6
    if(local_grid_x==0 && local_grid_y==Gy-1){
        vector<int> send_test_upper;
        vector<int> recv_test_upper(local_size_x+1);
        for (int i = 0; i < local_size_x+1; i++) {
            send_test_upper.push_back(local_data[i+(local_size_x+1)]);
        }
        MPI_Sendrecv(send_test_upper.data(), local_size_x+1, MPI_INT, rank-Gx, 0, recv_test_upper.data(), local_size_x+1, MPI_INT, rank-Gx, 0, MPI_COMM_WORLD, &st);
        for(int i=0; i<recv_test_upper.size(); i++){
            local_data[i]=recv_test_upper[i];
        }
    }
    if(local_grid_x<Gx-1 && local_grid_x>0 && local_grid_y==Gy-1){
        vector<int> send_test_upper;
        vector<int> recv_test_upper(local_size_x+2);
        for (int i = 0; i < local_size_x+2; i++) {
            send_test_upper.push_back(local_data[i+(local_size_x+2)]);
        }
        MPI_Sendrecv(send_test_upper.data(), local_size_x+2, MPI_INT, rank-Gx, 0, recv_test_upper.data(), local_size_x+2, MPI_INT, rank-Gx, 0, MPI_COMM_WORLD, &st);
        for(int i=0; i<recv_test_upper.size(); i++){
            local_data[i]=recv_test_upper[i];
        }
    }
    //case8
    if(local_grid_x==Gx-1 && local_grid_y==Gy-1){
        vector<int> send_test_upper;
        vector<int> recv_test_upper(local_size_x+1);
        for (int i = 0; i < local_size_x+1; i++) {
            send_test_upper.push_back(local_data[i+(local_size_x+1)]);
        }
        MPI_Sendrecv(send_test_upper.data(), local_size_x+1, MPI_INT, rank-Gx, 0, recv_test_upper.data(), local_size_x+1, MPI_INT, rank-Gx, 0, MPI_COMM_WORLD, &st);
        for(int i=0; i<recv_test_upper.size(); i++){
            local_data[i]=recv_test_upper[i];
        }
    }
}

void send_recv_data_y_double(vector<double> &local_data, int rank, int local_grid_x, int local_grid_y, int local_size_x, int local_size_y, int Gx, int Gy)
{
    //のりしろの通信
    MPI_Status st;
    ////case0
    if(local_grid_x==0 && local_grid_y==0){
        vector<double> send_test_lower;
        vector<double> recv_test_lower(local_size_x+1);
        for (int i = 0; i < local_size_x+1; i++) {
            send_test_lower.push_back(local_data[(local_size_x+1)*(local_size_y-1)+i]);
        }
        MPI_Sendrecv(send_test_lower.data(), local_size_x+1, MPI_DOUBLE, rank+Gx, 0, recv_test_lower.data(), local_size_x+1, MPI_DOUBLE, rank+Gx, 0, MPI_COMM_WORLD, &st);
        for(int i=0; i<recv_test_lower.size(); i++){
            local_data[(local_size_x+1)*local_size_y+i]=recv_test_lower[i];
        }
    }
    ////case1
    if(local_grid_y==0 && local_grid_x<Gx-1 && local_grid_x>0){
        vector<double> send_test_lower;
        vector<double> recv_test_lower(local_size_x+2);
        for (int i = 0; i < local_size_x+2; i++) {
            send_test_lower.push_back(local_data[(local_size_x+2)*(local_size_y-1)+i]);
        }
        MPI_Sendrecv(send_test_lower.data(), local_size_x+2, MPI_DOUBLE, rank+Gx, 0, recv_test_lower.data(), local_size_x+2, MPI_DOUBLE, rank+Gx, 0, MPI_COMM_WORLD, &st);
        for(int i=0; i<recv_test_lower.size(); i++){
            local_data[(local_size_x+2)*local_size_y+i]=recv_test_lower[i];
        }
    }
    ////case2
    if(local_grid_x==Gx-1 && local_grid_y==0){
        vector<double> send_test_lower;
        vector<double> recv_test_lower(local_size_x+1);
        for (int i = 0; i < local_size_x+1; i++) {
            send_test_lower.push_back(local_data[(local_size_x+1)*(local_size_y-1)+i]);
        }
        MPI_Sendrecv(send_test_lower.data(), local_size_x+1, MPI_DOUBLE, rank+Gx, 0, recv_test_lower.data(), local_size_x+1, MPI_DOUBLE, rank+Gx, 0, MPI_COMM_WORLD, &st);
        for(int i=0; i<recv_test_lower.size(); i++){
            local_data[(local_size_x+1)*local_size_y+i]=recv_test_lower[i];
        }
        
    }
    //case6
    if(local_grid_x==0 && local_grid_y==Gy-1){
        vector<double> send_test_upper;
        vector<double> recv_test_upper(local_size_x+1);
        for (int i = 0; i < local_size_x+1; i++) {
            send_test_upper.push_back(local_data[i+(local_size_x+1)]);
        }
        MPI_Sendrecv(send_test_upper.data(), local_size_x+1, MPI_DOUBLE, rank-Gx, 0, recv_test_upper.data(), local_size_x+1, MPI_DOUBLE, rank-Gx, 0, MPI_COMM_WORLD, &st);
        for(int i=0; i<recv_test_upper.size(); i++){
            local_data[i]=recv_test_upper[i];
        }
    }
    ////case7
    if(local_grid_x<Gx-1 && local_grid_x>0 && local_grid_y==Gy-1){
        vector<double> send_test_upper;
        vector<double> recv_test_upper(local_size_x+2);
        for (int i = 0; i < local_size_x+2; i++) {
            send_test_upper.push_back(local_data[i+(local_size_x+2)]);
        }
        MPI_Sendrecv(send_test_upper.data(), local_size_x+2, MPI_DOUBLE, rank-Gx, 0, recv_test_upper.data(), local_size_x+2, MPI_DOUBLE, rank-Gx, 0, MPI_COMM_WORLD, &st);
        for(int i=0; i<recv_test_upper.size(); i++){
            local_data[i]=recv_test_upper[i];
        }
    }
    //case8
    if(local_grid_x==Gx-1 && local_grid_y==Gy-1){
        vector<double> send_test_upper;
        vector<double> recv_test_upper(local_size_x+1);
        for (int i = 0; i < local_size_x+1; i++) {
            send_test_upper.push_back(local_data[i+(local_size_x+1)]);
        }
        MPI_Sendrecv(send_test_upper.data(), local_size_x+1, MPI_DOUBLE, rank-Gx, 0, recv_test_upper.data(), local_size_x+1, MPI_DOUBLE, rank-Gx, 0, MPI_COMM_WORLD, &st);
        for(int i=0; i<recv_test_upper.size(); i++){
            local_data[i]=recv_test_upper[i];
        }
    }
}

void calc_dxdr(int ic, std::vector<std::vector<double>> node, std::vector<std::vector<int>> element, std::vector<std::vector<double>> &dxdr, std::vector<std::vector<double>> dNdr)
{
  for(int k=0;k<2;k++){
    for(int l=0;l<2;l++){
      dxdr[k][l] = 0e0;
      for(int p=0;p<4;p++){
        dxdr[k][l] += dNdr[p][l] * node[element[ic][p]][k];
      }
    }
  }
}

void calc_inverse_matrix_2x2(std::vector<std::vector<double>> dxdr, std::vector<std::vector<double>> &drdx)
{
  double det = dxdr[0][0]*dxdr[1][1]-dxdr[0][1]*dxdr[1][0];
  drdx[0][0] = 1.0/det*dxdr[1][1];
  drdx[1][1] = 1.0/det*dxdr[0][0];
  drdx[0][1] = -1.0/det*dxdr[1][0];
  drdx[1][0] = -1.0/det*dxdr[0][1];
}

void calc_dNdx(vector<vector<double>> &dNdx, vector<vector<double>> dNdr, vector<vector<double>> drdx)
{
  for(int k=0; k<4; k++){
    for(int l=0; l<2; l++){
      dNdx[k][l] = 0.0;
      for(int p=0; p<2; p++){
        dNdx[k][l] += dNdr[k][p]*drdx[p][l];
      }
    }
  }
}

void remove_tmp_node(vector<double> &send_buf, vector<double> local_data, int local_grid_x, int local_grid_y, int local_size_x, int local_size_y, int Gx, int Gy)
{
    //case0
    if(local_grid_x==0 && local_grid_y==0){
        for(int i=0; i<local_size_y; i++){
            for(int j=0; j<local_size_x; j++){
                send_buf[i*local_size_x+j]=local_data[i*(local_size_x+1)+j];
            }
        }
    }
    //case1
    if(local_grid_y==0 && local_grid_x<Gx-1 && local_grid_x>0){
        for(int i=0; i<local_size_y; i++){
            for(int j=0; j<local_size_x; j++){
                send_buf[i*local_size_x+j]=local_data[i*(local_size_x+2)+j+1];
            }
        }
    }
    //case2
    else if(local_grid_x==Gx-1 && local_grid_y==0){
        for(int i=0; i<local_size_y; i++){
            for(int j=0; j<local_size_x; j++){
                send_buf[i*local_size_x+j]=local_data[i*(local_size_x+1)+j+1];
            }
        }
    }
    //case3
    else if(local_grid_x==0 && local_grid_y<Gy-1 && local_grid_y>0){
        for(int i=0; i<local_size_y; i++){
            for(int j=0; j<local_size_x; j++){
                send_buf[i*local_size_x+j]=local_data[(i+1)*(local_size_x+1)+j];
            }
        }
    }
    //case4
    else if(local_grid_x<Gx-1 && local_grid_x>0 && local_grid_y>0 && local_grid_y<Gy-1){
        for(int i=0; i<local_size_y; i++){
            for(int j=0; j<local_size_x; j++){
                send_buf[i*local_size_x+j]=local_data[(i+1)*(local_size_x+2)+j+1];
            }
        }
    }
    //case5
    else if(local_grid_x==Gx-1 && local_grid_y>0 && local_grid_y<Gy-1){
        for(int i=0; i<local_size_y; i++){
            for(int j=0; j<local_size_x; j++){
                send_buf[i*local_size_x+j]=local_data[(i+1)*(local_size_x+1)+j+1];
            }
        }
    }
    //case6
    else if(local_grid_x==0 && local_grid_y==Gy-1){
        for(int i=0; i<local_size_y; i++){
            for(int j=0; j<local_size_x; j++){
                send_buf[i*local_size_x+j]=local_data[(i+1)*(local_size_x+1)+j];
            }
        }
    }
    //case7
    else if(local_grid_x<Gx-1 && local_grid_x>0 && local_grid_y==Gy-1){
        for(int i=0; i<local_size_y; i++){
            for(int j=0; j<local_size_x; j++){
                send_buf[i*local_size_x+j]=local_data[(i+1)*(local_size_x+2)+j+1];
            }
        }
    }
    //case8
    else if(local_grid_x==Gx-1 && local_grid_y==Gy-1){
        for(int i=0; i<local_size_y; i++){
            for(int j=0; j<local_size_x; j++){
                send_buf[i*local_size_x+j]=local_data[(i+1)*(local_size_x+1)+j+1];
            }
        }
    }
}
void export_vtu(const std::string &file, vector<vector<int>> element, vector<vector<double>> node, vector<double> C)
{
    FILE *fp;
  if ((fp = fopen(file.c_str(), "w")) == NULL)
  {
    cout << file << " open error" << endl;
    exit(1);
  }

  fprintf(fp, "<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt32\">\n");
  fprintf(fp, "<UnstructuredGrid>\n");
  fprintf(fp, "<Piece NumberOfPoints= \"%d\" NumberOfCells= \"%d\" >\n", node.size(), element.size());
  fprintf(fp, "<Points>\n");
  int offset = 0;
  fprintf(fp, "<DataArray type=\"Float64\" Name=\"Position\" NumberOfComponents=\"3\" format=\"appended\" offset=\"%d\"/>\n",offset);
  offset += sizeof(int) + sizeof(double) * node.size() * 3;
  fprintf(fp, "</Points>\n");

  fprintf(fp, "<Cells>\n");
  fprintf(fp, "<DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">\n");
  for (int i = 0; i < element.size(); i++){
    for (int j = 0; j < element[i].size(); j++) fprintf(fp, "%d ", element[i][j]);
    fprintf(fp, "\n");
  }
  fprintf(fp, "</DataArray>\n");
  fprintf(fp, "<DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">\n");
  int num = 0;
  for (int i = 0; i < element.size(); i++)
  {
    num += element[i].size();
    fprintf(fp, "%d\n", num);
  }
  fprintf(fp, "</DataArray>\n");
  fprintf(fp, "<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n");
  for (int i = 0; i < element.size(); i++) fprintf(fp, "%d\n", 9);
    
  fprintf(fp, "</DataArray>\n");
  fprintf(fp, "</Cells>\n");

  fprintf(fp, "<PointData>\n");
  fprintf(fp, "<DataArray type=\"Float64\" Name=\"pressure[Pa]\" NumberOfComponents=\"1\" format=\"appended\" offset=\"%d\"/>\n",offset);
  offset += sizeof(int) + sizeof(double) * element.size();

  

  fprintf(fp, "</PointData>\n");

  fprintf(fp, "<CellData>\n");
  fprintf(fp, "</CellData>\n");
  fprintf(fp, "</Piece>\n");
  fprintf(fp, "</UnstructuredGrid>\n");
  fprintf(fp, "<AppendedData encoding=\"raw\">\n");
  fprintf(fp, "_");
  fclose(fp);

  ofstream ofs;
  ofs.open(file.c_str(), ios::out | ios::app | ios_base::binary);
  double *data_d = new double[node.size()*3];
  num = 0;
  int size=0;
  for (int ic = 0; ic < node.size(); ic++){
    data_d[num] = node[ic][0];
    num++;
    data_d[num] = node[ic][1];
    num++;
    data_d[num] = 0.0;
    num++;
  }
  size=sizeof(double)*node.size()*3;
  ofs.write((char *)&size, sizeof(size));
  ofs.write((char *)data_d, size);

  num=0;
  for (int ic = 0; ic < node.size(); ic++){
      data_d[num]   = C[ic];
      num++;
  }
  size=sizeof(double)*node.size();
  ofs.write((char *)&size, sizeof(size));
  ofs.write((char *)data_d, size);

  delete data_d;

  ofs.close();

  if ((fp = fopen(file.c_str(), "a")) == NULL)
  {
    cout << file << " open error" << endl;
    exit(1);
  }
  fprintf(fp, "\n</AppendedData>\n");
  fprintf(fp, "</VTKFile>\n");
  fclose(fp);
}
void change_order(vector<double> &last_data, vector<double> recvbuf, int local_size_x, int local_size_y, int point_num, int Gx)
{
    int local_box=local_size_x*local_size_y;
    for(int i=0; i<point_num; i++){
        for(int j=0; j<point_num; j++){
            int local_judge_x = j/local_size_x;
            int local_judge_y = i/local_size_y;
            int rank_tmp = Gx*local_judge_y+local_judge_x;
            last_data[i*point_num+j]=recvbuf[local_box*rank_tmp+local_size_x*(i%local_size_y)+(j%local_size_x)];
        }
    }   
}
int main(int argc, char **argv)
{
    clock_t start = clock();
    double dt = 0.01;
    int time_step = 36000;
    double diffusion_coefficient = 1e-2;
    vector<vector<double>> node;
    vector<vector<int>> element;
    int point_num = 64;
    set_field(node, element, point_num, 0.1);

    MPI_Init(&argc, &argv);
    int rank, procs;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &procs);
    int d2[2] = {};
    MPI_Dims_create(procs, 2, d2);

    int local_grid_x = rank%d2[0];
    int local_grid_y = rank/d2[0];
    int local_size_x = point_num/d2[0];
    int local_size_y = point_num/d2[1];
    vector<int> local_data_global;
    int corner=local_grid_y*(local_size_y*point_num)+local_size_x*(local_grid_x);
    set_global_node_number(local_data_global,local_grid_x,local_grid_y,local_size_x,local_size_y,corner,point_num,d2[0],d2[1]);    
    vector<double> local_data;
    int numofNode;
    //x方向通信
    send_recv_data_x_int(local_data_global,rank,local_grid_x,local_grid_y,local_size_x,local_size_y,d2[0],d2[1]);
    //y方向通信
    send_recv_data_y_int(local_data_global,rank,local_grid_x,local_grid_y,local_size_x,local_size_y,d2[0],d2[1]);
    vector<vector<int>> global_element;
    vector<vector<int>> local_element;
    
    //case0
    if(local_grid_x==0 && local_grid_y==0){
        local_data.resize((local_size_x+1)*(local_size_y+1));
        numofNode= (local_size_x+1)*(local_size_y+1);
        for(int i=0; i<local_size_y; i++){
            for(int j=0; j<local_size_x; j++){
                vector<int> global_element_tmp;
                global_element_tmp.push_back(local_data_global[(local_size_x+1)*i+j]);
                global_element_tmp.push_back(local_data_global[(local_size_x+1)*i+j]+local_size_x*d2[0]);
                global_element_tmp.push_back(local_data_global[(local_size_x+1)*i+j]+local_size_x*d2[0]+1);
                global_element_tmp.push_back(local_data_global[(local_size_x+1)*i+j]+1);
                global_element.push_back(global_element_tmp);
            }
        }
        for(int i=0; i<local_size_y; i++){
            for(int j=0; j<local_size_x; j++){
                vector<int> local_element_tmp;
                local_element_tmp.push_back((local_size_x+1)*i+j);
                local_element_tmp.push_back((local_size_x+1)*i+j+local_size_x+1);
                local_element_tmp.push_back((local_size_x+1)*i+j+local_size_x+2);
                local_element_tmp.push_back((local_size_x+1)*i+j+1);
                local_element.push_back(local_element_tmp);
            }
        }
    }
    //case1
    if(local_grid_y==0 && local_grid_x<d2[0]-1 && local_grid_x>0){
        local_data.resize((local_size_x+2)*(local_size_y+1));
        numofNode= (local_size_x+2)*(local_size_y+1);
        for(int i=0; i<local_size_y; i++){
            for(int j=0; j<local_size_x+1; j++){
                vector<int> global_element_tmp;
                global_element_tmp.push_back(local_data_global[(local_size_x+2)*i+j]);
                global_element_tmp.push_back(local_data_global[(local_size_x+2)*i+j]+(local_size_x)*d2[0]);
                global_element_tmp.push_back(local_data_global[(local_size_x+2)*i+j]+(local_size_x)*d2[0]+1);
                global_element_tmp.push_back(local_data_global[(local_size_x+2)*i+j]+1);
                global_element.push_back(global_element_tmp);
            }
        }
        for(int i=0; i<local_size_y; i++){
            for(int j=0; j<local_size_x+1; j++){
                vector<int> local_element_tmp;
                local_element_tmp.push_back((local_size_x+2)*i+j);
                local_element_tmp.push_back((local_size_x+2)*i+j+(local_size_x+1)+1);
                local_element_tmp.push_back((local_size_x+2)*i+j+(local_size_x+1)+2);
                local_element_tmp.push_back((local_size_x+2)*i+j+1);
                local_element.push_back(local_element_tmp);
            }
        }
    }
    //case2
    if(local_grid_x==d2[0]-1 && local_grid_y==0){
        local_data.resize((local_size_x+1)*(local_size_y+1));
        numofNode= (local_size_x+1)*(local_size_y+1);
        for(int i=0; i<local_size_y; i++){
            for(int j=0; j<local_size_x; j++){
                vector<int> global_element_tmp;
                global_element_tmp.push_back(local_data_global[(local_size_x+1)*i+j]);
                global_element_tmp.push_back(local_data_global[(local_size_x+1)*i+j]+local_size_x*d2[0]);
                global_element_tmp.push_back(local_data_global[(local_size_x+1)*i+j]+local_size_x*d2[0]+1);
                global_element_tmp.push_back(local_data_global[(local_size_x+1)*i+j]+1);
                global_element.push_back(global_element_tmp);
            }
        }
        for(int i=0; i<local_size_y; i++){
            for(int j=0; j<local_size_x; j++){
                vector<int> local_element_tmp;
                local_element_tmp.push_back((local_size_x+1)*i+j);
                local_element_tmp.push_back((local_size_x+1)*i+j+local_size_x+1);
                local_element_tmp.push_back((local_size_x+1)*i+j+local_size_x+2);
                local_element_tmp.push_back((local_size_x+1)*i+j+1);
                local_element.push_back(local_element_tmp);
            }
        }
    }
    //case6
    if(local_grid_x==0 && local_grid_y==d2[1]-1){
        local_data.resize((local_size_x+1)*(local_size_y+1));
        numofNode= (local_size_x+1)*(local_size_y+1);
        for(int i=0; i<local_size_y; i++){
            for(int j=0; j<local_size_x; j++){
                vector<int> global_element_tmp;
                global_element_tmp.push_back(local_data_global[(local_size_x+1)*i+j]);
                global_element_tmp.push_back(local_data_global[(local_size_x+1)*i+j]+local_size_x*d2[0]);
                global_element_tmp.push_back(local_data_global[(local_size_x+1)*i+j]+local_size_x*d2[0]+1);
                global_element_tmp.push_back(local_data_global[(local_size_x+1)*i+j]+1);
                global_element.push_back(global_element_tmp);
            }
        }
        for(int i=0; i<local_size_y; i++){
            for(int j=0; j<local_size_x; j++){
                vector<int> local_element_tmp;
                local_element_tmp.push_back((local_size_x+1)*i+j);
                local_element_tmp.push_back((local_size_x+1)*i+j+local_size_x+1);
                local_element_tmp.push_back((local_size_x+1)*i+j+local_size_x+2);
                local_element_tmp.push_back((local_size_x+1)*i+j+1);
                local_element.push_back(local_element_tmp);
            }
        }
    }
    //case7
    if(local_grid_x<d2[0]-1 && local_grid_x>0 && local_grid_y==d2[1]-1){
        local_data.resize((local_size_x+2)*(local_size_y+1));
        numofNode= (local_size_x+2)*(local_size_y+1);
        for(int i=0; i<local_size_y; i++){
            for(int j=0; j<local_size_x+1; j++){
                vector<int> global_element_tmp;
                global_element_tmp.push_back(local_data_global[(local_size_x+2)*i+j]);
                global_element_tmp.push_back(local_data_global[(local_size_x+2)*i+j]+(local_size_x)*d2[0]);
                global_element_tmp.push_back(local_data_global[(local_size_x+2)*i+j]+(local_size_x)*d2[0]+1);
                global_element_tmp.push_back(local_data_global[(local_size_x+2)*i+j]+1);
                global_element.push_back(global_element_tmp);
            }
        }
        for(int i=0; i<local_size_y; i++){
            for(int j=0; j<local_size_x+1; j++){
                vector<int> local_element_tmp;
                local_element_tmp.push_back((local_size_x+2)*i+j);
                local_element_tmp.push_back((local_size_x+2)*i+j+(local_size_x+1)+1);
                local_element_tmp.push_back((local_size_x+2)*i+j+(local_size_x+1)+2);
                local_element_tmp.push_back((local_size_x+2)*i+j+1);
                local_element.push_back(local_element_tmp);
            }
        }
    }
    //case8
    if(local_grid_x==d2[0]-1 && local_grid_y==d2[1]-1){
        local_data.resize((local_size_x+1)*(local_size_y+1));
        numofNode= (local_size_x+1)*(local_size_y+1);
        for(int i=0; i<local_size_y; i++){
            for(int j=0; j<local_size_x; j++){
                vector<int> global_element_tmp;
                global_element_tmp.push_back(local_data_global[(local_size_x+1)*i+j]);
                global_element_tmp.push_back(local_data_global[(local_size_x+1)*i+j]+local_size_x*d2[0]);
                global_element_tmp.push_back(local_data_global[(local_size_x+1)*i+j]+local_size_x*d2[0]+1);
                global_element_tmp.push_back(local_data_global[(local_size_x+1)*i+j]+1);
                global_element.push_back(global_element_tmp);
            }
        }
        for(int i=0; i<local_size_y; i++){
            for(int j=0; j<local_size_x; j++){
                vector<int> local_element_tmp;
                local_element_tmp.push_back((local_size_x+1)*i+j);
                local_element_tmp.push_back((local_size_x+1)*i+j+local_size_x+1);
                local_element_tmp.push_back((local_size_x+1)*i+j+local_size_x+2);
                local_element_tmp.push_back((local_size_x+1)*i+j+1);
                local_element.push_back(local_element_tmp);
            }
        }
    }

    if(rank==0){
        local_data[0] = 1.0;
    }
    
    vector<vector<double>> gauss(4, vector<double>(2));
    gauss[0][0] = -0.577350296189626; gauss[0][1] = -0.577350296189626;
    gauss[1][0] = -0.577350296189626; gauss[1][1] = 0.577350296189626; 
    gauss[2][0] = 0.577350296189626; gauss[2][1] = -0.577350296189626; 
    gauss[3][0] = 0.577350296189626; gauss[3][1] = 0.577350296189626;

    vector<vector<double>> D;
    vector<vector<double>> mass;
    vector<double> mass_centralization;
    D.resize(numofNode);
    for(int i=0; i<numofNode; i++){
        D[i].resize(numofNode);
    }
    mass.resize(numofNode);
    for(int i=0; i<numofNode; i++){
        mass[i].resize(numofNode);
    }
    mass_centralization.resize(numofNode);
    
    for(int i=0; i<local_element.size(); i++){
        vector<double> N(4); 
        vector<vector<double>> dNdr(4, vector<double>(2));
        double volume=0.0;
        vector<vector<double>> element_D(4, vector<double>(4,0.0));
        vector<vector<double>> element_mass(4, vector<double>(4,0.0));
        vector<double> element_mass_centralization(4,0.0);
        for(int j=0; j<4; j++){
            ShapeFunction2D::C2D4_N(N,gauss[j][0],gauss[j][1]);
            ShapeFunction2D::C2D4_dNdr(dNdr,gauss[j][0],gauss[j][1]);
            vector<vector<double>> dxdr(2, vector<double>(2)), drdx(2, vector<double>(2));
            
            calc_dxdr(i, node, global_element, dxdr, dNdr);
            calc_inverse_matrix_2x2(dxdr, drdx);
            double detJ = -(dxdr[0][0] * dxdr[1][1]  - dxdr[1][0] * dxdr[0][1]);
            vector<vector<double>> dNdx(4, vector<double>(2, 0.0));
            calc_dNdx(dNdx, dNdr, drdx);
            for(int k=0; k<4; k++){
              for(int l=0; l<4; l++){
                for(int p=0; p<2; p++){
                  element_D[k][l] += dNdx[k][p]*dNdx[l][p];
                }
              }
            }
            for(int k=0; k<4; k++){
              for(int l=0; l<4; l++){
                element_mass[k][l] = N[k] * N[l];
              }
            }
            for(int k=0; k<4; k++){
              for(int l=0; l<4; l++){
                D[local_element[i][k]][local_element[i][l]] += diffusion_coefficient * element_D[k][l] * detJ;
                mass[local_element[i][k]][local_element[i][l]] += element_mass[k][l] * detJ;
              }
            }
        }
    }
    for(int i=0; i<mass.size(); i++){
      mass_centralization[i]=0.0;
      for(int j=0; j<mass[i].size(); j++){
        mass_centralization[i] += mass[i][j];
      }
    }
    for(int i=0; i<time_step; i++){
        if(rank==0) cout << i << endl;
        vector<double> DC(numofNode,0.0);
        vector<double> DcR(numofNode,0.0);
        vector<double> MDcR(numofNode,0.0);
        vector<double> MDC(numofNode,0.0);

        for(int j=0; j<numofNode; j++){
            for(int k=0; k<numofNode; k++){
                DC[j] += D[j][k] * (local_data[k]);
            }
        }

        for(int j=0; j<numofNode; j++){
            MDC[j] = 1.0/mass_centralization[j]*DC[j];
            local_data[j] = local_data[j] - dt * MDC[j];
        }
        if(rank==0){
            local_data[0] = 1.0;
        }

        send_recv_data_x_double(local_data, rank, local_grid_x, local_grid_y, local_size_x, local_size_y, d2[0], d2[1]);
        send_recv_data_y_double(local_data,rank,local_grid_x,local_grid_y,local_size_x,local_size_y,d2[0],d2[1]);
        
        if(i%1000==0){
            vector<double> send_buf(local_size_x*local_size_y);
            remove_tmp_node(send_buf,local_data,local_grid_x,local_grid_y,local_size_x,local_size_y,d2[0],d2[1]);

            vector<double> recvbuf;
            if (rank == 0) {
                recvbuf.resize(local_size_x * local_size_y * procs);
            }
            MPI_Gather(send_buf.data(), local_size_x * local_size_y, MPI_DOUBLE, recvbuf.data(), local_size_x * local_size_y, MPI_DOUBLE, 0, MPI_COMM_WORLD);
            if(rank==0){
                vector<double> last_data(local_size_x*local_size_y*procs);
                change_order(last_data,recvbuf,local_size_x,local_size_y,point_num,d2[0]);
                int local_box=local_size_x*local_size_y;
                string filename = "test_mpi/test_"+to_string(i/1000)+".vtu";
                export_vtu(filename,element,node,last_data);
            }
        }
    }
    clock_t end = clock();   
    std::cout << "duration = " << (double)(end - start) / CLOCKS_PER_SEC << "sec.\n";
    MPI_Finalize();
}