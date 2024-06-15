// サポート構造を構築する

#include<iostream>
#include<random>
#include<vector>
#include<Eigen/Core>
#include<Eigen/Geometry>
#include<Eigen/Dense>
#include<cmath>
#include "readobj.hpp"
#include "overhang.hpp"

using namespace std;
using namespace Eigen;

// 底面(x,z平面)でのグリッドセルを作る(入力：モデルの座標,グリッドの分割数)
vector<Vector3d> grid_cell(vector<Vector3d> V,int n){
    vector<Vector3d> gc;
    
    double x_min=100,x_max=-100,z_min=100,z_max=-100,y_min=100;
    for(int i=0;i<V.size();i++){
        if(x_min>V[i](0)){
            x_min=V[i](0);
        }
        if(z_min>V[i](2)){
           z_min=V[i](2); 
        }
        if(x_max<V[i](0)){
            x_max=V[i](0);
        }
        if(z_max<V[i](2)){
            z_max=V[i](2);
        }
        if(y_min>V[i](1)){
            y_min=V[i](1);
        }
    }

    Vector3d g;
    // cout<<"x_min="<<x_min<<" y_min="<<y_min<<endl;

    double x_length=x_max-x_min; //横の長さ
    double y_length=z_max-z_min; //縦の長さ
    for(double i=z_min;i<=z_max;i+=y_length/n){
        for(double j=x_min;j<=x_max;j+=x_length/n){
            // cout<<"("<<i<<" ,"<<j<<")"<<endl;
            g(0)=j;
            g(1)=y_min;
            g(2)=i;
            gc.push_back(g);
        }
    }

    // for(int i=0;i<V.size();i++){
    //     cout<<gc[i]<<endl;
    // }

    return gc;
}

// 点だけobjファイル出力
void obj_out(vector<Vector3d> vert,const char* file_name){
    FILE* f;
    f=fopen(file_name,"w");
    
    for(int i=0;i<vert.size();i++){
        cout<<vert[i](0)<<","<<vert[i](1)<<endl;
        fprintf(f,"v %lf %lf %lf\n",vert[i](0),vert[i](1),vert[i](2));
    }
}


bool point_in_out(Vector3d oh,Vector2d s){

}

// サポート構築(Lattice)
vector<Vector3d> L_support(Vector3d s){
    
}

// サポート構築(Tree)
vector<Vector3d> T_support(){

}

int main(int argc,char* argv[])
{   
    string fn;
    FILE *fp;
    vector<Vector3d> V;
    vector<Vector3i> F;
    ro reado;

    if(argc==1){
        cout<<"入力してください\n";
        exit(1);
    }else{
        fn=argv[1];
        fp=fopen(argv[1],"r");
    }
    reado.readPoint(V,F,fn);
    cout<<"読み込み完了\n";
    // for(int i=0;i<V.size();i++){
    //     cout<<V[i]<<endl;
    // }
    // for(int i=0;i<F.size();i++){
    //     cout<<F[i]<<endl;
    // }

    // x,y平面グリッドセルの作成
    vector<Vector3d> gc=grid_cell(V,10);
    cout<<"グリッドセルの作成完了\n";
    obj_out(gc,"grid_cell.obj");

    // // 面のオーバハングを調べる
    // cout<<"面のオーバーハングを調べる"<<F.size()<<endl;
    // double angle=0.6108652381980153; //閾値:cura35度(0.6108652381980153)，45度：0.7853981633974483
    // Vector3d direction(0,-1,0); //造形方向ベクトル
    // vector<int> num; //オーバーハングかどうか(1 or 0)
    // for(int i=0;i<F.size();i++){
    //     if(overh(angle,direction,V[F[i](0)],V[F[i](1)],V[F[i](2)])) num.push_back(1);
    //     else num.push_back(0);
    // }

    // // vtk出力，結果の可視化
    // FILE* fpv;
    // fpv=fopen("overhang.vtk","w");
    // fprintf(fpv,"# vtk DataFile Version 2.0\n");
    // fprintf(fpv,"Title of my super cool VTK file\nASCII\nDATASET UNSTRUCTURED_GRID\n");
    // fprintf(fpv,"POINTS %d float\n",(int)V.size());
    // for(int i=0;i<V.size();i++){
    //     fprintf(fpv,"%lf %lf %lf\n",V[i](0),V[i](1),V[i](2));
    // }
    // fprintf(fpv,"CELLS %d %d\n",(int)F.size(),(int)F.size()*4);
    // for(int i=0;i<F.size();i++){
    //     fprintf(fpv,"3 %d %d %d\n",F[i](0),F[i](1),F[i](2));
    // }
    // fprintf(fpv,"CELL_TYPES %d\n",(int)F.size());
    // for(int i=0;i<F.size();i++){
    //     fprintf(fpv,"5\n");
    // }
    // fprintf(fpv,"CELL_DATA %d\n",(int)F.size());
    // fprintf(fpv,"SCALARS my_face_attr int 1\nLOOKUP_TABLE default\n");
    // for(int i=0;i<num.size();i++){
    //     fprintf(fpv,"%d\n",num[i]);
    // }
    // cout<<"numsize="<<num.size()<<endl;
}

// 面のクラス分け
// オーバハングの面だけ残す
// 隣の面がオーバハングか探索する
// https://docs.google.com/spreadsheets/d/1N3K9jsz0g74o5LToadMPtwWfvS1vS7ZjCbrNMR23uYA/edit?usp=sharing