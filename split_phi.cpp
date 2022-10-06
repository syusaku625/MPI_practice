#include<fstream>
#include<iostream>
#include<string>
#include<vector>

using namespace std;

void export_vtu(const std::string &file, vector<double> output_value, vector<vector<int>> element, vector<vector<double>> node)
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
  fprintf(fp, "</PointData>\n");

  fprintf(fp, "<CellData>\n");
  fprintf(fp, "<DataArray type=\"Float64\" Name=\"pressure[Pa]\" NumberOfComponents=\"1\" format=\"appended\" offset=\"%d\"/>\n",offset);
  offset += sizeof(int) + sizeof(double) * element.size();
  fprintf(fp, "</CellData>\n");
  fprintf(fp, "</Piece>\n");
  fprintf(fp, "</UnstructuredGrid>\n");
  fprintf(fp, "<AppendedData encoding=\"raw\">\n");
  fprintf(fp, "_");
  fclose(fp);

  cout << "test" << endl;

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
  delete data_d;
  cout << "test" << endl;

  double *data_d2 = new double[element.size()];
  num=0;
  for (int ic = 0; ic < element.size(); ic++){
    data_d2[num] = output_value[ic];
    cout << output_value[ic] << endl;
    num++;
  }
  cout << "test" << endl;
  size=sizeof(double)*element.size();
  ofs.write((char *)&size, sizeof(size));
  ofs.write((char *)data_d2, size);
  delete data_d2;
  cout << "test2" << endl;

  ofs.close();

  if ((fp = fopen(file.c_str(), "a")) == NULL)
  {
    cout << file << " open error" << endl;
    exit(1);
  }
  fprintf(fp, "\n</AppendedData>\n");
  fprintf(fp, "</VTKFile>\n");
  fclose(fp);
  cout << "test" << endl;
}

void set_field(vector<vector<double>> &node, vector<vector<int>> &element, int nx, int ny, double dx)
{
    node.resize(nx*ny);
    for(int i=0; i<node.size(); i++){
        node[i].resize(2);
    }
    for(int i=0; i<ny; i++){
        for(int j=0; j<nx; j++){
            node[i*nx+j][0]= j*dx;
            node[i*nx+j][1]= i*dx;
        }
    }
    element.resize((nx-1)*(ny-1));
    for(int i=0; i<ny-1; i++){
        for(int j=0; j<nx-1; j++){
            element[i*(nx-1)+j].push_back(i*nx+j);
            element[i*(nx-1)+j].push_back(i*nx+j+1);
            element[i*(nx-1)+j].push_back(i*nx+j+(nx+1));
            element[i*(nx-1)+j].push_back(i*nx+j+nx);   
        }
    }
}

int main()
{
    int numOfElm = 127*127;
    vector<double> phi;
    ifstream ifs("fluid_phi.dat");
    string str,tmp;
    phi.resize(numOfElm);
    for(int i=0; i<numOfElm; i++){
      getline(ifs,str);
      phi[i] = stod(str);
    }
    ifs.close();

    vector<vector<double>> phi_trans(127, vector<double>(127));
    for(int i=0; i<127; i++){
        for(int j=0; j<127; j++){
            phi_trans[i][j]=phi[i*127+j];
        }
    }

    vector<double> phi1(32*64);
    vector<double> phi2(32*64);
    vector<double> phi3(32*64);
    vector<double> phi4(32*64);
    vector<double> phi5(32*64);
    vector<double> phi6(32*64);
    vector<double> phi7(32*64);
    vector<double> phi8(32*64);

    for(int i=0; i<64; i++){
        for(int j=0; j<32; j++){
            phi1[i*63+j]=phi_trans[i][j];
        }
    }
    for(int i=0; i<64; i++){
        for(int j=31; j<63; j++){
            phi2[i*63+j]=phi_trans[i][j];
        }
    }
    //for(int i=0; i<64; i++){
    //    for(int j=62; j<95; j++){
    //        phi3[i*64+j]=phi_trans[i][j];
    //    }
    //}
    //for(int i=0; i<64; i++){
    //    for(int j=94; j<127; j++){
    //        phi4[i*64+j]=phi_trans[i][j];
    //    }
    //}

    vector<vector<int>> element;
    vector<vector<double>> node;
    cout << "test" << endl;
    set_field(node,element,33,65,0.1);
    cout << element.size() << " " << phi1.size() << endl;
    cout << "test" << endl;
    string file="phi_test.vtu";
    export_vtu(file,phi2,element,node);   
    cout << "end" << endl;
}