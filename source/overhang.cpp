// オーバハングを検出するプログラム

#include<iostream>
#include<Eigen/Core>
#include<Eigen/Geometry>
#include<cmath>
#include "readobj.hpp"
#include "overhang.hpp"

using namespace std;
using namespace Eigen;

// 法線ベクトルを求める(引数：三点の座標)
Vector3d normal(Vector3d o,Vector3d a,Vector3d b){
    Vector3d A=a-o;
    Vector3d B=b-o;
    return A.cross(B);
}

// オーバハングかどうか調べる(引数：角度の閾値，造形する方向下向きベクトル，選んだ面の三点)オーバハングの場合trueを返す
bool overh(double angle,Vector3d direction,Vector3d o,Vector3d a,Vector3d b){
    Vector3d n=normal(o,a,b);
    // cout<<"n="<<n<<endl;
    Vector3d A=a-o;
    Vector3d B=b-o;
    double shita=acos(direction.dot(n)/(direction.norm()*n.norm())); //造形方向ベクトルと法線ベクトルのなす角
    // cout<<shita<<endl;
    if(angle>shita){
        // cout<<"true\n";
        return true;
    }else{
        // cout<<"false\n";
        return false;
    }
}

// int main(int argc,char* argv[])
// {   
//     string fn;
//     FILE *fp;
//     vector<Vector3d> V;
//     vector<Vector3i> F;
//     ro reado;

//     if(argc==1){
//         cout<<"入力してください\n";
//         exit(1);
//     }else{
//         fn=argv[1];
//         fp=fopen(argv[1],"r");
//     }
//     reado.readPoint(V,F,fn);
//     // for(int i=0;i<V.size();i++){
//     //     cout<<V[i]<<endl;
//     // }
//     // for(int i=0;i<F.size();i++){
//     //     cout<<F[i]<<endl;
//     // }

//     // 面のオーバハングを調べる
//     double angle=0.5235987755982988; //閾値:cura35度(0.6108652381980153)，45度：0.7853981633974483
//     Vector3d direction(0,-1,0); //造形方向ベクトル
//     vector<int> num; //オーバーハングかどうか(1 or 0)
//     for(int i=0;i<F.size();i++){
//         if(overh(angle,direction,V[F[i](0)],V[F[i](1)],V[F[i](2)])) num.push_back(1);
//         else num.push_back(0);
//     }

//     // vtk出力，結果の可視化
//     FILE* fpv;
//     fpv=fopen("overhang.vtk","w");
//     fprintf(fpv,"# vtk DataFile Version 2.0\n");
//     fprintf(fpv,"Title of my super cool VTK file\nASCII\nDATASET UNSTRUCTURED_GRID\n");
//     fprintf(fpv,"POINTS %d float\n",(int)V.size());
//     for(int i=0;i<V.size();i++){
//         fprintf(fpv,"%lf %lf %lf\n",V[i](0),V[i](1),V[i](2));
//     }
//     fprintf(fpv,"CELLS %d %d\n",(int)F.size(),(int)F.size()*4);
//     for(int i=0;i<F.size();i++){
//         fprintf(fpv,"3 %d %d %d\n",F[i](0),F[i](1),F[i](2));
//     }
//     fprintf(fpv,"CELL_TYPES %d\n",(int)F.size());
//     for(int i=0;i<F.size();i++){
//         fprintf(fpv,"5\n");
//     }
//     fprintf(fpv,"CELL_DATA %d\n",(int)F.size());
//     fprintf(fpv,"SCALARS my_face_attr int 1\nLOOKUP_TABLE default\n");
//     for(int i=0;i<num.size();i++){
//         fprintf(fpv,"%d\n",num[i]);
//     }
//     cout<<"numsize="<<num.size()<<endl;
// }