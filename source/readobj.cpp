// objを読み込むプログラム

#include<iostream>
#include<Eigen/Core>
#include<Eigen/Geometry>
#include<fstream>
#include<sstream>
#include<vector>
#include<string>
#include "readobj.hpp"

using namespace std;
using namespace Eigen;

// 座標と面を読み込む()
void ro::readPoint(vector<Vector3d> &vertex,vector<Vector3i> &face,string fn){
    cout<<"読み込み中\n";
    string str,dm;
    ifstream ifs(fn);

    Vector3d v;
    Vector3i f;

    while(getline(ifs,str)){
        stringstream ss(str);
        if(ss.str()[0]=='v'){
            ss>>dm;
            ss>>v(0)>>v(1)>>v(2);
            vertex.push_back(v);
        }
        if(ss.str()[0]=='f'){
            ss>>dm;
            ss>>f(0)>>f(1)>>f(2);
            Vector3i fa={1,1,1};
            face.push_back(f-fa);
        }
    }
}

// int main(int argc,char* argv[])
// {
//     string fn;
//     FILE *fp;
//     if(argc==1){
//         cout<<"入力してください\n";
//         exit(1);
//     }else{
//         fn=argv[1];
//         fp=fopen(argv[1],"r");
//     }

//     vector<Vector3d> V;
//     vector<Vector3i> F;
//     readPoint(V,F,fn);

//     // for(int i=0;i<V.size();i++){
//     //     cout<<V[i]<<endl;
//     // }
//     // for(int i=0;i<F.size();i++){
//     //     cout<<F[i]<<endl;
//     // }
// }