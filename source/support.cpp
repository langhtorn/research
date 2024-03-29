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

//六角形のクラス
struct hexagon
{
    int index; //六角形の番号
    vector<int> hi; //今見ている六角形の回りの六角形のindex
};

// 60°の計算
Vector3d make_etriangle(Vector3d p,Vector3d pd,Vector3d n,double length){
    Vector3d d;
    MatrixXd A(3,3);
    A.row(0)=p-pd;
    A.row(1)=pd-p;
    A.row(2)=n;
    Vector3d B(3);
    B(0)=(1/2)*length*length+(pd.dot(p))-pd.dot(pd);
    B(1)=(1/2)*length*length+(pd.dot(p))-p.dot(p);
    B(2)=0;
    d=A.colPivHouseholderQr().solve(B);
    return d;
}

// 点の内外判定
bool point_in_out(Vector3d oh,Vector2d s){

}

// ハニカムサンプリング(入力：最小包含の座標大小,一辺の長さ)
vector<Vector2d> honeycomb(Vector2d min,Vector2d max,double r){
    vector<Vector2d> sample;
    double i=min(0),j=min(1);
    while(j<max(1)){ //最大のy座標
        i=min(0);
        while(i<max(0)){ //最大のx座標
            Vector2d center{i,j};
            vector<Vector2d> mh=make_hexagon(center,r);
            sample.insert(sample.end(),mh.begin(),mh.end());
            i+=r*3;
        }
        j+=r*sin(1.04719755);
    }
    return sample;
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
    // for(int i=0;i<V.size();i++){
    //     cout<<V[i]<<endl;
    // }
    // for(int i=0;i<F.size();i++){
    //     cout<<F[i]<<endl;
    // }

    // 面のオーバハングを調べる
    double angle=0.5235987755982988; //閾値:cura35度(0.6108652381980153)，45度：0.7853981633974483
    Vector3d direction(0,-1,0); //造形方向ベクトル
    vector<int> num; //オーバーハングかどうか(1 or 0)
    for(int i=0;i<F.size();i++){
        if(overh(angle,direction,V[F[i](0)],V[F[i](1)],V[F[i](2)])) num.push_back(1);
        else num.push_back(0);
    }

    // vtk出力，結果の可視化
    FILE* fpv;
    fpv=fopen("overhang.vtk","w");
    fprintf(fpv,"# vtk DataFile Version 2.0\n");
    fprintf(fpv,"Title of my super cool VTK file\nASCII\nDATASET UNSTRUCTURED_GRID\n");
    fprintf(fpv,"POINTS %d float\n",(int)V.size());
    for(int i=0;i<V.size();i++){
        fprintf(fpv,"%lf %lf %lf\n",V[i](0),V[i](1),V[i](2));
    }
    fprintf(fpv,"CELLS %d %d\n",(int)F.size(),(int)F.size()*4);
    for(int i=0;i<F.size();i++){
        fprintf(fpv,"3 %d %d %d\n",F[i](0),F[i](1),F[i](2));
    }
    fprintf(fpv,"CELL_TYPES %d\n",(int)F.size());
    for(int i=0;i<F.size();i++){
        fprintf(fpv,"5\n");
    }
    fprintf(fpv,"CELL_DATA %d\n",(int)F.size());
    fprintf(fpv,"SCALARS my_face_attr int 1\nLOOKUP_TABLE default\n");
    for(int i=0;i<num.size();i++){
        fprintf(fpv,"%d\n",num[i]);
    }
    // cout<<"numsize="<<num.size()<<endl;
}

// 面のクラス分け
// オーバハングの面だけ残す
// 隣の面がオーバハングか探索する