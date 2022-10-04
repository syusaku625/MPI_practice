#include<vector>
#include<iostream>
#include"shapefunction.hpp"
#include<fstream>
#include <time.h>

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

void calc_dxdr(int ic, std::vector<std::vector<double>> node, std::vector<std::vector<int>> element, std::vector<std::vector<double>> &dxdr, std::vector<std::vector<double>> dNdr)
{
  int k,l,p;
  for(k=0;k<2;k++){
    for(l=0;l<2;l++){
      dxdr[k][l] = 0e0;
      for(p=0;p<4;p++){
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

  fstream ofs;
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

int main(int argc, char **argv)
{
    clock_t start = clock();
    double dt = 0.01;
    int time_step = 36000;
    double diffusion_coefficient = 1e-2;
    vector<vector<double>> node;
    vector<vector<int>> element;
    int point_num = 32;
    set_field(node, element, point_num, 0.1);
    vector<vector<double>> gauss(4, vector<double>(2));
    gauss[0][0] = -0.577350296189626; gauss[0][1] = -0.577350296189626;
    gauss[1][0] = -0.577350296189626; gauss[1][1] = 0.577350296189626; 
    gauss[2][0] = 0.577350296189626; gauss[2][1] = -0.577350296189626; 
    gauss[3][0] = 0.577350296189626; gauss[3][1] = 0.577350296189626;

    vector<vector<double>> D;
    vector<vector<double>> mass;
    vector<double> mass_centralization;
    vector<double> C(node.size());
    D.resize(node.size());
    for(int i=0; i<node.size(); i++){
        D[i].resize(node.size());
    }
    mass.resize(node.size());
    for(int i=0; i<node.size(); i++){
        mass[i].resize(node.size());
    }
    mass_centralization.resize(node.size());

    for(int i=0; i<element.size(); i++){
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
            calc_dxdr(i, node, element, dxdr, dNdr);
            calc_inverse_matrix_2x2(dxdr, drdx);
            double detJ = dxdr[0][0] * dxdr[1][1]  - dxdr[1][0] * dxdr[0][1];
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
                D[element[i][k]][element[i][l]] += diffusion_coefficient * element_D[k][l] * detJ;
                mass[element[i][k]][element[i][l]] += element_mass[k][l] * detJ;
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
        cout << i << endl;
        vector<double> DC(node.size(),0.0);
        vector<double> DcR(node.size(),0.0);
        vector<double> MDcR(node.size(),0.0);
        vector<double> MDC(node.size(),0.0);
    
        for(int i=0; i<node.size(); i++){
            for(int j=0; j<node.size(); j++){
                DC[i] += D[i][j] * (C[j]);
            }
        }

        for(int i=0; i<node.size(); i++){
            MDC[i] = 1.0/mass_centralization[i]*DC[i];
            C[i] = C[i] - dt * MDC[i];
        }

        C[0]=1.0;
        if(i%1000==0){
            string filename = "test/test_"+to_string(i/1000)+".vtu";
            cout << filename << endl;
            export_vtu(filename,element,node,C);
        }
    }
    clock_t end = clock();   
    std::cout << "duration = " << (double)(end - start) / CLOCKS_PER_SEC << "sec.\n";
}
