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

struct Model{
    vector<Vector3d> V;
    vector<Vector3i> F;
};

vector<int> NoRemovableFaces; // 除去不能な面のインデックス

void ExportVTK(const Model& G, const vector<int>& NoRemovableFaces, const string& filename) {
    ofstream file(filename);
    if (!file) return;

    // VTK ヘッダー
    file << "# vtk DataFile Version 3.0\n";
    file << "VTK output\n";
    file << "ASCII\n";
    file << "DATASET POLYDATA\n";

    // 頂点の出力
    file << "POINTS " << G.V.size() << " float\n";
    for (const auto& v : G.V) {
        file << v.x() << " " << v.y() << " " << v.z() << "\n";
    }

    // 面の出力
    file << "POLYGONS " << G.F.size() << " " << G.F.size() * 4 << "\n";
    for (const auto& f : G.F) {
        file << "3 " << f.x() << " " << f.y() << " " << f.z() << "\n";
    }

    // 面の色（スカラー値として出力）
    file << "CELL_DATA " << G.F.size() << "\n";
    file << "SCALARS FaceColors float 1\n";
    file << "LOOKUP_TABLE default\n";

    for (size_t i = 0; i < G.F.size(); ++i) {
        file << (find(NoRemovableFaces.begin(), NoRemovableFaces.end(), i) != NoRemovableFaces.end() ? 1.0f : 0.0f) << "\n";
    }

    file.close();
}


// 球面サンプリング
// // 辺の中点を計算する関数
// Vector3d judge::findMidpoint(Vector3d A,Vector3d B){
//     Vector3d midpoint=(A+B)/2.0;
//     return midpoint.normalized()*r;
// }

// // 正二十面体の初期解
// std::vector<Eigen::Vector3i> generate_icosahedron_faces() {
//     return {
//         { 0, 11, 5 }, { 0, 5, 1 }, { 0, 1, 7 }, { 0, 7, 10 }, { 0, 10, 11 },
//         { 1, 5, 9 }, { 5, 11, 4 }, { 11, 10, 2 }, { 10, 7, 6 }, { 7, 1, 8 },
//         { 3, 9, 4 }, { 3, 4, 2 }, { 3, 2, 6 }, { 3, 6, 8 }, { 3, 8, 9 },
//         { 4, 9, 5 }, { 2, 4, 11 }, { 6, 2, 10 }, { 8, 6, 7 }, { 9, 8, 1 }
//     };
// }

// std::vector<Eigen::Vector3d> generate_icosahedron_vertices(const Eigen::Vector3d& center_point) {
//     std::vector<Eigen::Vector3d> vertices;
//     double phi = (1.0 + std::sqrt(5.0)) / 2.0;

//     std::vector<Eigen::Vector3d> unit_vertices = {
//         { -1,  phi,  0 }, {  1,  phi,  0 }, { -1, -phi,  0 }, {  1, -phi,  0 },
//         {  0, -1,  phi }, {  0,  1,  phi }, {  0, -1, -phi }, {  0,  1, -phi },
//         {  phi,  0, -1 }, {  phi,  0,  1 }, { -phi,  0, -1 }, { -phi,  0,  1 }
//     };

//     for (const auto& v : unit_vertices) {
//         vertices.push_back(center_point + (v.normalized() * r));
//     }

//     return vertices;
// }


// // サブディビジョンする関数
// void judge::subdivide(int n,Vector3d centerpoint){
//     samplingV=generate_icosahedron_vertices();
//     samplingF=generate_icosahedron_faces();
//     vector<Vector3d> newV;
//     vector<Vector3i> newF;
//     for(int j=0;j<n;j++){
//         for(int i=0;i<samplingF.size();i++){
//             int v1=F[i](0);
//             int v2=F[i](1);
//             int v3=F[i](2);

//             // 中点の計算
//             Vector3d m12v=findMidpoint(V[v1],V[v2]);
//             Vector3d m23v=findMindpoint(V[v2],V[v3]);
//             Vector3d m31v=findMindpoint(V[v3],V[v1]);
//             int m12=newV.size(); // 新しい頂点インデックス
//             newV.push_back(m12v);
//             int m23=newV.size();
//             newV.push_back(m23v);
//             int m31=newV.size();
//             newV.push_back(m31v);

//             // 新しい三角形を生成
//             newF.push_back({v1,m12,m31});
//             newF.push_back({v2,m23,m12});
//             newF.push_back({v3,m31,m23});
//             newF.push_back({m12,m23,m31});
//         }
//         F=newF;
//         V=newV;
//     }
// }

int inum=0;

// OBJファイル形式で3Dメッシュを出力する関数
    void visualizeMeshToObj(const vector<Vector3d>& vertices, const vector<Vector3i>& faces, const string& filename) {
        ofstream objFile(filename);
        
        if (!objFile.is_open()) {
            cerr << "Failed to open file for writing: " << filename << endl;
            return;
        }

        // 頂点データの書き込み (v x y z)
        for (const auto& vertex : vertices) {
            objFile << "v " << vertex.x() << " " << vertex.y() << " " << vertex.z() << endl;
        }

        // 面データの書き込み (f v1 v2 v3)
        for (const auto& face : faces) {
            objFile << "f " << face.x() + 1 << " " << face.y() + 1 << " " << face.z() + 1 << endl;
        }

        objFile.close();
        cout << "OBJ file written to: " << filename << endl;
    }


// 三角錐の中に円錐が完全に含まれるか
// 三角錐の作成
Model pyramid(vector<Vector3d> triangle,vector<Vector3d> up){

    Model pyramid;

    pyramid.V=triangle;
    Vector3d g=(up[0]+up[1]+up[2])/3.0;
    pyramid.V.push_back(g);

    pyramid.F.push_back(Vector3i(0,1,3));
    pyramid.F.push_back(Vector3i(1,2,3));
    pyramid.F.push_back(Vector3i(2,0,3));
    pyramid.F.push_back(Vector3i(1,0,2));

    return pyramid;
}
// 四面体の体積を計算する関数
double computeTetrahedronVolume(const Vector3d& A, const Vector3d& B, const Vector3d& C, const Vector3d& D) {
    return std::abs((B - A).dot((C - A).cross(D - A))) / 6.0;
}
// 三角形の面積を求める関数
double computeTriangleArea(const Vector3d& A, const Vector3d& B, const Vector3d& C) {
    return 0.5 * ((B - A).cross(C - A)).norm();
}
// 三角錐の内接球の重心座標
Vector3d computeCentroid(Model p){
    Vector3d A=p.V[3];
    Vector3d B=p.V[0];
    Vector3d C=p.V[1];
    Vector3d D=p.V[2];
    double Sa,Sb,Sc,Sd;
    Sa=computeTriangleArea(B,C,D);
    Sb=computeTriangleArea(A,C,D);
    Sc=computeTriangleArea(A,D,B);
    Sd=computeTriangleArea(A,B,C);

    return (A * Sa + B * Sb + C * Sc + D * Sd) / (Sa + Sb + Sc + Sd);

}
// 内接球の半径
double computeTriangleInradius(Model p){
    Vector3d A=p.V[3];
    Vector3d B=p.V[0];
    Vector3d C=p.V[1];
    Vector3d D=p.V[2];
    double Sa,Sb,Sc,Sd;
    Sa=computeTriangleArea(B,C,D);
    Sb=computeTriangleArea(A,C,D);
    Sc=computeTriangleArea(A,D,B);
    Sd=computeTriangleArea(A,B,C);

    double V=computeTetrahedronVolume(A,B,C,D);
    return V*(3.0/(Sa+Sb+Sc+Sd));
}
// 内接球の円錐の角度を求める
double ConeAngle(double radius,double l){
    return asin(radius/l);
}
// 球面三角形と元モデルの面からサポートの除去可能判定をする
bool isSupportRemoval(vector<Vector3d> Si,vector<Vector3d> MG){
    // 三角錐の作成
    Model p=pyramid(Si,MG);
    if(inum==396) visualizeMeshToObj(p.V,p.F,"sannkakusui.obj");

    // 内接球を求める
    Vector3d centroid=computeCentroid(p); // 中心座標
    double radius=computeTriangleInradius(p); //半径

    double l=(p.V[3]-centroid).norm();
    double cangle=ConeAngle(radius,l);

    double toolangle=0.183261; // 10.5°
    if(cangle>=toolangle){
        return true;
    }else{
        return false;
    }

}

// アクセス可能性行列の読み込み
vector<vector<bool>> readAccessStatusFromFile(const string& filename) {
    ifstream file(filename);
    if (!file) {
        cerr << "Error: Could not open file " << filename << endl;
        return {};
    }

    vector<vector<bool>> AS;
    string line;

    while (getline(file, line)) {
        istringstream iss(line);
        vector<bool> row;
        int value;
        
        while (iss >> value) {
            row.push_back(value);  // 0 または 1 を bool に変換
        }

        AS.push_back(row);
    }

    file.close();
    return AS;
}

// 各オーバーハング面に対して除去可能かみて，不可能になった面の面積の総計を求める
double RemovalSupportArea(vector<vector<bool>> AS,Model S,Model G){
    double area=0;
    for(int i=0;i<G.F.size();i++){
        bool i_flg=false;
        for(int j=0;j<AS.size();j++){
            if(AS[j][i]){
                vector<Vector3d> Si={S.V[S.F[j][0]],S.V[S.F[j][1]],S.V[S.F[j][2]]};
                vector<Vector3d> Gi={G.V[G.F[i][0]],G.V[G.F[i][1]],G.V[G.F[i][2]]};
                if(isSupportRemoval(Si,Gi)){
                    cout<<"除去可能\n";
                    i_flg=true;
                    break;
                }
            }
        }
        if(i_flg==false){
            area+=computeTriangleArea(G.V[G.F[i][0]],G.V[G.F[i][1]],G.V[G.F[i][2]]);
            NoRemovableFaces.push_back(i);
        }
        inum++;
    }

    return area;
}


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
    vector<int> overhung_index;
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

    // アクセス可能性行列の読み込み
    vector<vector<bool>> AS=readAccessStatusFromFile("accesstatusinfo.txt");

    Model S;
    reado.readPoint(S.V,S.F,"sphericaltriangle.obj");
    Model G=Model{V,F};
    cout<<"モデル作成\n";

    double sp_area=RemovalSupportArea(AS,S,G);

    cout<<"sp_area="<<sp_area<<endl;

    ExportVTK(G,NoRemovableFaces,"noremovalface.vtk");
    
}