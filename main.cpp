#include<mpi.h>
#include<vector>
#include<iostream>
using namespace std;

void set_field(vector<vector<double>> &node, vector<vector<int>> &element, int n, double dx)
{
    node.resize((n+1)*(n+1));
    for(int i=0; i<node.size(); i++){
        node[i].resize(2);
    }
    for(int i=0; i<n+1; i++){
        for(int j=0; j<n+1; j++){
            node[i*(n+1)+j][0]= j*dx;
            node[i*(n+1)+j][1]= i*dx;
        }
    }
    element.resize(n*n);
    for(int i=0; i<n; i++){
        for(int j=0; j<n; j++){
            element[i*n+j].push_back(i*(n+1)+j);
            element[i*n+j].push_back(i*(n+1)+j+1);
            element[i*n+j].push_back(i*(n+1)+j+(n+2));
            element[i*n+j].push_back(i*(n+1)+j+(n+1));   
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
                else if(j==local_size_x){
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

void remove_tmp_node(vector<int> &send_buf, vector<int> local_data, int local_grid_x, int local_grid_y, int local_size_x, int local_size_y, int Gx, int Gy)
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

void change_order(vector<int> &last_data, vector<int> recvbuf, int local_size_x, int local_size_y, int point_num, int Gx)
{
    int local_box=local_size_x*local_size_y;
    for(int i=0; i<point_num; i++){
        for(int j=0; j<point_num; j++){
            int local_judge_x = j/local_size_x;
            int local_judge_y = i/local_size_y;
            int rank_tmp = Gx*local_judge_y+local_judge_x;
            last_data[i*point_num+j]=recvbuf[local_box*rank_tmp+local_size_y*(i%local_size_y)+(j%local_size_x)];
        }
    }   
}

void send_recv_data_x(vector<int> &local_data, int rank, int local_grid_x, int local_grid_y, int local_size_x, int local_size_y, int Gx, int Gy)
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
    //case2
    if(local_grid_x==Gx-1 && local_grid_y==0){
        vector<int> send_test_left;
        vector<int> recv_test_left(local_size_y);
        for (int i = 0; i < local_size_y; i++) {
            send_test_left.push_back(local_data[(local_size_x+1)*i+1]);
        }
        MPI_Sendrecv(send_test_left.data(), local_size_y, MPI_INT, rank-1, 0, recv_test_left.data(), local_size_x, MPI_INT, rank-1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
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
        MPI_Sendrecv(send_test_right.data(), local_size_y, MPI_INT, rank+1, 0, recv_test_right.data(), local_size_y, MPI_INT, rank+1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        for(int i=0; i<recv_test_right.size(); i++){
            local_data[(local_size_x+1)*(i+1)+local_size_x]=recv_test_right[i];
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

void send_recv_data_y(vector<int> &local_data, int rank, int local_grid_x, int local_grid_y, int local_size_x, int local_size_y, int Gx, int Gy)
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
    //case2
    if(local_grid_x==Gx-1 && local_grid_y==0){
        vector<int> send_test_lower;
        vector<int> recv_test_lower(local_size_x+1);
        for (int i = 0; i < local_size_x+1; i++) {
            send_test_lower.push_back(local_data[(local_size_x+1)*(local_size_y-1)+i]);
        }
        MPI_Sendrecv(send_test_lower.data(), local_size_x+1, MPI_INT, rank+Gx, 0, recv_test_lower.data(), local_size_x+1, MPI_INT, rank+Gx, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
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
        MPI_Sendrecv(send_test_upper.data(), local_size_x+1, MPI_INT, rank-Gx, 0, recv_test_upper.data(), local_size_x+1, MPI_INT, rank-Gx, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
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

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);
    int rank, procs;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &procs);
    int d2[2] = {};
    MPI_Dims_create(procs, 2, d2);
    vector<vector<double>> node;
    vector<vector<int>> element;
    int point_num = 8;
    set_field(node, element, point_num, 0.1);
    int local_grid_x = rank%d2[0];
    int local_grid_y = rank/d2[0];
    int local_size_x = point_num/d2[0];
    int local_size_y = point_num/d2[1];
    vector<int> local_data;
    int corner=local_grid_y*(local_size_y*point_num)+local_size_x*(local_grid_x);
    set_global_node_number(local_data,local_grid_x,local_grid_y,local_size_x,local_size_y,corner,point_num,d2[0],d2[1]);    
    //x方向通信
    send_recv_data_x(local_data,rank,local_grid_x,local_grid_y,local_size_x,local_size_y,d2[0],d2[1]);
    //y方向通信
    send_recv_data_y(local_data,rank,local_grid_x,local_grid_y,local_size_x,local_size_y,d2[0],d2[1]);
    //のりしろの削除
    vector<int> send_buf(local_size_x*local_size_y);
    remove_tmp_node(send_buf,local_data,local_grid_x,local_grid_y,local_size_x,local_size_y,d2[0],d2[1]);

    vector<int> recvbuf;
    if (rank == 0) {
        recvbuf.resize(local_size_x * local_size_y * procs);
    }
    
    MPI_Gather(send_buf.data(), local_size_x * local_size_y, MPI_INT, recvbuf.data(), local_size_x * local_size_y, MPI_INT, 0, MPI_COMM_WORLD);
    if(rank==0){
        vector<int> last_data(local_size_x*local_size_y*procs);
        change_order(last_data,recvbuf,local_size_x,local_size_y,point_num,d2[0]);
        int local_box=local_size_x*local_size_y;
        for(int i=0; i<last_data.size(); i++){
            cout << last_data[i] << endl;
        }
    }
    MPI_Finalize();
}