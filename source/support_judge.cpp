//サポートの良否判定

#include<iostream>
#include<random>
#include<vector>
#include<Eigen/Core>
#include<Eigen/Geometry>
#include<Eigen/Dense>
#include<cmath>
#include<igl/copyleft/cgal/mesh_boolean.h>
#include "readobj.hpp"
#include "overhang.hpp"
#include "vtomat.hpp"
#include<igl/writeOBJ.h>
#include "support.hpp"
#include "support_judge.hpp"
#include<igl/copyleft/cgal/intersect_other.h>

using namespace std;
using namespace Eigen;

// 球面サンプリング
// 辺の中点を計算する関数
Vector3d judge::findMidpoint(Vector3d A,Vector3d B){
    Vector3d midpoint=(A+B)/2.0;
    return midpoint.normalized()*r;
}

// 正二十面体の初期解
std::vector<Eigen::Vector3i> generate_icosahedron_faces() {
    return {
        { 0, 11, 5 }, { 0, 5, 1 }, { 0, 1, 7 }, { 0, 7, 10 }, { 0, 10, 11 },
        { 1, 5, 9 }, { 5, 11, 4 }, { 11, 10, 2 }, { 10, 7, 6 }, { 7, 1, 8 },
        { 3, 9, 4 }, { 3, 4, 2 }, { 3, 2, 6 }, { 3, 6, 8 }, { 3, 8, 9 },
        { 4, 9, 5 }, { 2, 4, 11 }, { 6, 2, 10 }, { 8, 6, 7 }, { 9, 8, 1 }
    };
}

std::vector<Eigen::Vector3d> generate_icosahedron_vertices(const Eigen::Vector3d& center_point) {
    std::vector<Eigen::Vector3d> vertices;
    double phi = (1.0 + std::sqrt(5.0)) / 2.0;

    std::vector<Eigen::Vector3d> unit_vertices = {
        { -1,  phi,  0 }, {  1,  phi,  0 }, { -1, -phi,  0 }, {  1, -phi,  0 },
        {  0, -1,  phi }, {  0,  1,  phi }, {  0, -1, -phi }, {  0,  1, -phi },
        {  phi,  0, -1 }, {  phi,  0,  1 }, { -phi,  0, -1 }, { -phi,  0,  1 }
    };

    for (const auto& v : unit_vertices) {
        vertices.push_back(center_point + (v.normalized() * r));
    }

    return vertices;
}


// サブディビジョンする関数
void judge::subdivide(int n,Vector3d centerpoint){
    samplingV=generate_icosahedron_vertices();
    samplingF=generate_icosahedron_faces();
    vector<Vector3d> newV;
    vector<Vector3i> newF;
    for(int j=0;j<n;j++){
        for(int i=0;i<samplingF.size();i++){
            int v1=F[i](0);
            int v2=F[i](1);
            int v3=F[i](2);

            // 中点の計算
            Vector3d m12v=findMidpoint(V[v1],V[v2]);
            Vector3d m23v=findMindpoint(V[v2],V[v3]);
            Vector3d m31v=findMindpoint(V[v3],V[v1]);
            int m12=newV.size(); // 新しい頂点インデックス
            newV.push_back(m12v);
            int m23=newV.size();
            newV.push_back(m23v);
            int m31=newV.size();
            newV.push_back(m31v);

            // 新しい三角形を生成
            newF.push_back({v1,m12,m31});
            newF.push_back({v2,m23,m12});
            newF.push_back({v3,m31,m23});
            newF.push_back({m12,m23,m31});
        }
        F=newF;
        V=newV;
    }
}


// AccessStatus行列のfalseの割合が閾値を超えたら除去不能
vector<bool> determine

// // モデルAABBツリーと三角形メッシュとの交差判定
// bool judge::intersect_triangle(MatrixXd TV,MatrixXi TF){

//     MatrixXi IF; // 交差する三角形のペアが格納される
//     bool kari=true;

//     bool intersects=igl::copyleft::cgal::intersect_other(sp.VG_Mver,sp.FG_Mver,TV,TF,kari,IF);

//     if(intersects) return true;
//     else return false;

// }

// // 中心点から周囲の１点を求める(中心点と正規化された方向ベクトル)
// MatrixXd judge::four_point(Vector3d centerpoint,Vector3d direction){

//     MatrixXd fp(1,3);
//     Vector3d v0=direction.unitOrthogonal();
//     Vector3d v1=direction.cross(v0);
//     vector<Vector3d> p;

//     // 周囲の4点求める
//     double scale=1; //中心点からの距離
//     p.push_back(centerpoint+scale*v0+scale*v1);
//     p.push_back(centerpoint+scale*v0-scale*v1);
//     p.push_back(centerpoint-scale*v0+scale*v1);
//     p.push_back(centerpoint-scale*v0-scale*v1);

//     fp.row(0)=p[0].transpose();

//     for(int i=1;i<4;i++){
//         fp.conservativeResize(fp.rows()+1,fp.cols());
//         fp.row(fp.rows()-1)=p[i].transpose();
//     }
    
//     return fp;
// }

// // 削除できるサポートを探す
// void judge::delete_suppport(){
//     for(int i=0;i<sp.ohvn.size();i++){ //オーバーハング点の番号，サポート番号

//         for(int j=0;j<sp.rays_s[i].size();j++){ //始点を決める
//             // 道具ボクセルの作成
            
//             Vector3d direction(0,-1,0); // 初期方向ベクトルは造形方向下向きで
//             Vector3d centerpoint=sp.rays_s[i][j]; //始点を平面の中心点とする

//             direction.normalize(); //正規化して単位ベクトルにする
            
//             MatrixXd toolpoints; //ツールの頂点座標
//             MatrixXi toolface; //ツールの面座標

//             MatrixXd tp1=four_point(centerpoint,direction); //周囲の4点を求める

//             double L=1; //ツールの長さ
//             centerpoint=centerpoint+L*direction; // ツールボクセルの反対側

//             MatrixXd tp2=four_point(centerpoint,-direction); //下の周囲の4点
//             toolpoints(tp1.rows()+tp2.rows(),tp1.cols());
//             toolpoints<<tp1,
//                         tp2;

//             bool intersects=intersect_triangle(toolpoints,toolface); //交差判定 trueが交差してる
//         }
//     }
// }

int main(int argc,char* argv[])
{   
    string fn;
    FILE *fp;
    vector<Vector3d> V;
    vector<Vector3i> F;
    ro reado;
    spt sp;

    if(argc==1){
        cout<<"入力してください\n";
        exit(1);
    }else{
        fn=argv[1];
        fp=fopen(argv[1],"r");
    }
    reado.readPoint(V,F,fn);
    sp.VG=V;
    sp.FG=F;
    cout<<"読み込み完了\n";
    // for(int i=0;i<V.size();i++){
    //     cout<<V[i]<<endl;
    // }
    // for(int i=0;i<F.size();i++){
    //     cout<<F[i]<<endl;
    // }

    // x,y平面グリッドセルの作成
    sp.N=30; //グリッドの分割数
    sp.aabb(); //最小包含
    sp.grid_cell();
    cout<<"グリッドセルの作成完了\n";
    // obj_out(gc,"grid_cell.obj"); //obj出力

    double y_max=-100;
    for(int i=0;i<V.size();i++){ //y座標最大
        if(y_max<V[i](1)){
            y_max=V[i](1);
        }
    }

    // 面のオーバハングを調べる
    cout<<"面のオーバーハングを調べる"<<sp.FG.size()<<endl;
    double angle=0.5235987755983; //閾値:cura35度(0.6108652381980153)，45度：0.7853981633974483,30度0.5235987755983
    Vector3d direction(0,-1,0); //造形方向ベクトル
    vector<int> num; //その面がオーバーハングかどうか(1 or 0)
    vectir<int> overhung_index;
    for(int i=0;i<F.size();i++){
        if(overh(angle,direction,V[F[i](0)],V[F[i](1)],V[F[i](2)])){ 
            num.push_back(1);
            overhung_index.push_back(i);
        }else{
            num.push_back(0);
        } 
    }

    // オーバハング点の判定
    sp.oh_Fnum(angle,direction); // オーバーハング面の集合を作る
    // Vの射影
    for(int i=0;i<V.size();i++){
        Vector3d vss={V[i](0),sp.four[0](1),V[i](2)};
        sp.Vs.push_back(vss);
    }

    for(int i=0;i<sp.oh_fn.size();i++){
        Vector3i ohf={F[sp.oh_fn[i]](0),F[sp.oh_fn[i]](1),F[sp.oh_fn[i]](2)};
        sp.oh_F.push_back(ohf);
    }
    double h=y_max-sp.four[0](1);
    vector<int> oh_point=sp.oh_Vnum(sp.gc_V,h); //オーバーハング点の点番号
    
    // デバッグ
    cout<<"gc="<<sp.gc_V.size()<<" ohp="<<oh_point.size()<<endl;
    for(int i=0;i<oh_point.size();i++){
        Vector3d op=sp.gc_V[oh_point[i]];
        sp.ohp.push_back(op);
    }
    cout<<"ohpsize="<<sp.ohp.size()<<endl;

    sp.obj_out(sp.ohp,"oh_point.obj");

    // vector<Vector3d> sv; //サポート点情報
    // vector<Vector3i> sf; //サポート面情報
    // sp.l_height();
    // cout<<"h="<<h<<endl;
    // sp.L_support();
    // cout<<"sv="<<sv.size()<<" sf="<<sf.size()<<endl;
    // sp.obj_outf(sv,sf,"support_L.obj");

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