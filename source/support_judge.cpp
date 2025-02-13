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
#include<igl/ray_mesh_intersect.h>
#include<igl/per_face_normals.h>

using namespace std;
using namespace Eigen;

struct Model{
    vector<Vector3d> V;
    vector<Vector3i> F;
};

// VTKファイルを書き出す関数
void writeVTK(const string& filename, const MatrixXd& V, const MatrixXi& F, const vector<pair<int,int>>& oh_pair) {
    ofstream file(filename);
    if (!file) {
        cerr << "エラー: VTKファイルを開けません: " << filename << endl;
        return;
    }

    file << "# vtk DataFile Version 3.0\n";
    file << "VTK output from support_judge\n";
    file << "ASCII\n";
    file << "DATASET POLYDATA\n";

    // 頂点の書き出し
    file << "POINTS " << V.rows() << " float\n";
    for (int i = 0; i < V.rows(); i++) {
        file << V(i, 0) << " " << V(i, 1) << " " << V(i, 2) << "\n";
    }

    // 面の書き出し
    int num_faces = 0;
    for (const auto& p : oh_pair) {
        num_faces++;  // オーバーハング面
        if (p.second != -1) num_faces++;  // 相方の面
    }
    file << "POLYGONS " << num_faces << " " << num_faces * 4 << "\n";

    for (const auto& p : oh_pair) {
        file << "3 " << F(p.first, 0) << " " << F(p.first, 1) << " " << F(p.first, 2) << "\n";
        if (p.second != -1) {
            file << "3 " << F(p.second, 0) << " " << F(p.second, 1) << " " << F(p.second, 2) << "\n";
        }
    }

    // 色データの書き出し
    file << "CELL_DATA " << num_faces << "\n";
    file << "SCALARS FaceType int 1\n";
    file << "LOOKUP_TABLE default\n";

    for (const auto& p : oh_pair) {
        file << "1\n"; // オーバーハング面（赤）
        if (p.second != -1) {
            file << "2\n"; // 相方の面（青）
        }
    }

    file.close();
    cout << "✅ VTKファイルを書き出しました: " << filename << endl;
}

// VTKファイルを書き出す関数
void writeVTKoh(const string& filename, const MatrixXd& V, const MatrixXi& F, const vector<int>& oh_fn) {
    ofstream file(filename);
    if (!file) {
        cerr << "エラー: VTKファイルを開けません: " << filename << endl;
        return;
    }

    file << "# vtk DataFile Version 3.0\n";
    file << "VTK output for overhang faces\n";
    file << "ASCII\n";
    file << "DATASET POLYDATA\n";

    // 頂点の書き出し
    file << "POINTS " << V.rows() << " float\n";
    for (int i = 0; i < V.rows(); i++) {
        file << V(i, 0) << " " << V(i, 1) << " " << V(i, 2) << "\n";
    }

    // 面の書き出し
    file << "POLYGONS " << F.rows() << " " << F.rows() * 4 << "\n";
    for (int i = 0; i < F.rows(); i++) {
        file << "3 " << F(i, 0) << " " << F(i, 1) << " " << F(i, 2) << "\n";
    }

    // 色付けのための CELL_DATA を作成
    unordered_set<int> oh_faces(oh_fn.begin(), oh_fn.end()); // 高速検索用
    file << "CELL_DATA " << F.rows() << "\n";
    file << "SCALARS OverhangFace int 1\n";
    file << "LOOKUP_TABLE default\n";
    
    for (int i = 0; i < F.rows(); i++) {
        if (oh_faces.count(i)) {
            file << "1\n";  // オーバーハング面 → 赤色
        } else {
            file << "0\n";  // 通常の面 → 灰色
        }
    }

    file.close();
    cout << "✅ VTKファイルを書き出しました: " << filename << endl;
}

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
// 辺の中点を計算する関数
Vector3d findMidpoint(Vector3d A,Vector3d B){
    Vector3d midpoint=(A+B)/2.0;
    return midpoint.normalized();
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
        vertices.push_back(center_point + (v.normalized() ));
    }

    return vertices;
}


// サブディビジョンする関数
vector<Vector3d> subdivide(int n){
    Vector3d center_point={0,0,0};
    vector<Vector3d> V=generate_icosahedron_vertices(center_point);
    vector<Vector3i> F=generate_icosahedron_faces();
    vector<Vector3d> newV;
    vector<Vector3i> newF;
    for(int j=0;j<n;j++){
        for(int i=0;i<F.size();i++){
            int v1=F[i](0);
            int v2=F[i](1);
            int v3=F[i](2);

            // 中点の計算
            Vector3d m12v=findMidpoint(V[v1],V[v2]);
            Vector3d m23v=findMidpoint(V[v2],V[v3]);
            Vector3d m31v=findMidpoint(V[v3],V[v1]);
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
    return V;
}

int inum=0;
int jnum=0;

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
// 三角錐の頂点選択
Vector3d selectApex(vector<Vector3d> MG){
    Vector3d apex=MG[0];
    for(int i=1;i<3;i++){
        if(MG[i].y()>apex.y()){
            apex=MG[i];
        }
    }
    return apex;
}
// 三角錐の作成
Model pyramid(vector<Vector3d> triangle,vector<Vector3d> up){

    Model pyramid;

    pyramid.V=triangle;
    Vector3d g=selectApex(up);
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
// 正規化座標
double bari(double Sv,double Sa,double Sb,double Sc,double Sd){
    return Sv/(Sa+Sb+Sc+Sd);
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

    double rama,ramb,ramc,ramd;
    rama=bari(Sa,Sa,Sb,Sc,Sd);
    ramb=bari(Sb,Sa,Sb,Sc,Sd);
    ramc=bari(Sc,Sa,Sb,Sc,Sd);
    ramd=bari(Sd,Sa,Sb,Sc,Sd);

    return rama*A+ramb*B+ramc*C+ramd*D;
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
    // if(inum==287 && jnum==0) visualizeMeshToObj(p.V,p.F,"sannkakusui.obj");

    // 内接球を求める
    Vector3d centroid=computeCentroid(p); // 中心座標
    double radius=computeTriangleInradius(p); //半径

    double l=(p.V[3]-centroid).norm();
    double cangle=ConeAngle(radius,l);
    // if(inum==287 && jnum==0) cout<<"cangle="<<cangle<<endl;

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
double RemovalSupportArea(vector<vector<bool>> AS,Model S,Model G,vector<pair<int,int>> oh,vector<int>& NoRemovableFaces){
    NoRemovableFaces.clear();
    double area=0;
    for(int i=0;i<oh.size();i++){
        
        bool i_flg=false;
        bool i_flg_aikata=false;
        int ohnum=oh[i].first;
        int aikata=oh[i].second;
        // cout<<"情報取得\n";
        for(int j=0;j<AS.size();j++){

            // オーバーハング面の除去判定
            if(AS[j][ohnum]){
                if(ohnum==514) cout<<"中\n";
                // cout<<"ohnum="<<ohnum<<endl;
                vector<Vector3d> Si={S.V[S.F[j][0]],S.V[S.F[j][1]],S.V[S.F[j][2]]};
                vector<Vector3d> Gi={G.V[G.F[ohnum][0]],G.V[G.F[ohnum][1]],G.V[G.F[ohnum][2]]};
                
                if(isSupportRemoval(Si,Gi)){
                    // cout<<i<<" 除去可能\n";

                    i_flg=true;
                }
                jnum++;
            }
            // 相方の除去判定
            if(aikata==-1 && i_flg==true){
                i_flg_aikata=true;
                break;
            }
            // オーバーハング面が除去可能で相方が地面でないとき
            if(i_flg && aikata!=-1 && AS[j][aikata]){
                // cout<<"aikata="<<aikata<<endl;
                vector<Vector3d> Si={S.V[S.F[j][0]],S.V[S.F[j][1]],S.V[S.F[j][2]]};
                vector<Vector3d> Gi={G.V[G.F[aikata][0]],G.V[G.F[aikata][1]],G.V[G.F[aikata][2]]};
                
                if(isSupportRemoval(Si,Gi)){
                    // cout<<i<<" 除去可能\n";

                    i_flg_aikata=true;
                    break;
                }
                jnum++;
            }

        }
        jnum=0;
        // cout<<"除去判定おわり\n";
        // 除去不可能になった面積の集計
        if(i_flg==false || i_flg_aikata==false){
            area+=computeTriangleArea(G.V[G.F[ohnum][0]],G.V[G.F[ohnum][1]],G.V[G.F[ohnum][2]]);
             NoRemovableFaces.push_back(ohnum);

            if(aikata!=-1){
                area+=computeTriangleArea(G.V[G.F[aikata][0]],G.V[G.F[aikata][1]],G.V[G.F[aikata][2]]);
                NoRemovableFaces.push_back(aikata);
            }
           
        }
        
        // cout<<"面積の集計\n";
        inum++;
    }
    // cout<<"NoRemovalFaces="<<NoRemovableFaces.size()<<endl;
    

    return area;
}

// オーバーハング面と向かいあっている面のペアを作る
pair<int,int> findfirstIntersection(Vector3d oh_point,int fnum,Vector3d direction,MatrixXd mv,MatrixXi mf,MatrixXd mf_normals){
    vector<igl::Hit> hits;

    direction.normalize();
direction += Vector3d(1e-9, 1e-9, 1e-9); // 誤差補正


    // レイとメッシュの交差を計算
    igl::ray_mesh_intersect(oh_point,direction,mv,mf,hits);

    if(hits.empty()) return {fnum,-1};

    // 最も近い交差点の面を取得
    int closestFace=hits[0].id;
    double min_t=hits[0].t;

    for(const auto& hit : hits){

        if(hit.id==fnum) continue;

        // 交差した面の法線を取得
        Vector3d faceNormal = mf_normals.row(hit.id).normalized();

        // レイの方向と法線の向きを確認
        if (faceNormal.dot(direction) >= 0) continue;  // 逆向きなら無視

        if(hit.t<min_t){
            min_t=hit.t;
            closestFace=hit.id;
        }
    }

    return {fnum,closestFace};

}

// VTKファイルを書き出す関数（1つのベクトル用）
void writeVTKVector(const string& filename, const Vector3d& start, const Vector3d& direction) {
    ofstream file(filename);
    if (!file) {
        cerr << "エラー: VTKファイルを開けません: " << filename << endl;
        return;
    }

    file << "# vtk DataFile Version 3.0\n";
    file << "VTK Vector Output\n";
    file << "ASCII\n";
    file << "DATASET POLYDATA\n";

    // 1つの点（ベクトルの始点）を出力
    file << "POINTS 1 float\n";
    file << start(0) << " " << start(1) << " " << start(2) << "\n";

    // ベクトルデータの出力
    file << "POINT_DATA 1\n";
    file << "VECTORS direction float\n";
    file << direction(0) << " " << direction(1) << " " << direction(2) << "\n";

    file.close();
    cout << "✅ VTKファイルを書き出しました: " << filename << endl;
}

// ラジアンを度に変換
double radToDeg(double rad) {
    return rad * (180.0 / M_PI);
}

int main(int argc,char* argv[])
{   
    string fn;
    FILE *fp;
    vector<Vector3d> V;
    vector<Vector3i> F;
    ro reado;
    spt sp;
    vt vtom;

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

    // アクセス可能性行列の読み込み
    vector<vector<bool>> AS=readAccessStatusFromFile("accesstatusinfo.txt");

    Model S;
    reado.readPoint(S.V,S.F,"sphericaltriangle.obj");
    Model G=Model{V,F};
    MatrixXd mv;
    MatrixXi mf;
    vtom.vmat(G.V,G.F,mv,mf);
    cout<<"モデル作成\n";

    vector<Vector3d> buildDirections=subdivide(1);
    sort(buildDirections.begin(), buildDirections.end(), [](const Vector3d& a, const Vector3d& b) {
    return a.norm() < b.norm();
    });

    vector<int> NoremovalFaces_min;
    vector<int> NoremovalFaces_max;

    // 面のオーバハングを調べる
    cout<<"面のオーバーハングを調べる"<<sp.FG.size()<<endl;
    double angle=0.5235987755983; //閾値:cura35度(0.6108652381980153)，45度：0.7853981633974483,30度0.5235987755983
    double sp_area_min=10000; // 最小除去不能面積
    Vector3d min_vector;
    vector<int> Noremoval_min;
    double sp_area_max=-1; // 最大除去不能面積
    Vector3d max_vector;
    vector<int> Noremoval_max;
    vector<pair<int,int>> max_overhung;
    cout<<"builddirection="<<buildDirections.size()<<endl;
    
    for(int i=0;i<buildDirections.size();i++){
        // cout<<"造形方向回転開始\n";
        cout<<i<<endl;
        sp.ohp.clear();
        sp.oh_F.clear();
        sp.oh_fn.clear();
        Vector3d direction=buildDirections[i]; // 造形方向ベクトル
        direction.y()=direction.y()*(-1);
        // Vector3d direction(0,-1,0);
        vector<int> num; //その面がオーバーハングかどうか(1 or 0)
        vector<int> overhung_index;
        for(int j=0;j<F.size();j++){
            if(overh(angle,direction,V[F[j](0)],V[F[j](1)],V[F[j](2)])){ 
                num.push_back(1);
                overhung_index.push_back(j);
            }else{
                num.push_back(0);
            } 
        }
        // cout<<"オーバーハング面の作成\n";

        // オーバハング点の判定
        sp.oh_Fnum(angle,direction); // オーバーハング面の集合を作る
        // cout<<"オーバーハング面の面番号集合作成\n";
        // // Vの射影
        // for(int j=0;j<V.size();j++){
        //     Vector3d vss={V[j](0),sp.four[0](1),V[j](2)};
        //     sp.Vs.push_back(vss);
        // }

        for(int j=0;j<sp.oh_fn.size();j++){
            Vector3i ohf={F[sp.oh_fn[j]](0),F[sp.oh_fn[j]](1),F[sp.oh_fn[j]](2)};
            sp.oh_F.push_back(ohf);
        }

        // オーバーハング面と向かいあっている面のペアを作る
        MatrixXd mf_normals;
        igl::per_face_normals(mv,mf,mf_normals);
        vector<pair<int,int>> oh_pair;
        for(int j=0;j<sp.oh_fn.size();j++){
            Vector3i f=sp.FG[sp.oh_fn[j]]; // オーバハング面の頂点インデックス
            Vector3d centroid=(sp.VG[f(0)]+sp.VG[f(1)]+sp.VG[f(2)])/3.0;
            pair<int,int> oh_p=findfirstIntersection(centroid,sp.oh_fn[j],direction,mv,mf,mf_normals);
            oh_pair.push_back(oh_p);
        }
        
        vector<int> NoRemovableFaces;
        double sp_area=RemovalSupportArea(AS,S,G,oh_pair,NoRemovableFaces);
        // cout<<"除去可能な場所の判定\n";

        if(sp_area<sp_area_min){
            sp_area_min=sp_area;
            direction.y()=-direction.y();
            min_vector=direction;
            NoremovalFaces_min=NoRemovableFaces;
            max_overhung=oh_pair;
        }

        if(sp_area>sp_area_max){
            sp_area_max=sp_area;
            direction.y()=-direction.y();
            max_vector=direction;
            NoremovalFaces_max=NoRemovableFaces;
            // max_overhung=oh_pair;
        }

        // cout<<"noremovalface="<<NoRemovableFaces.size();

        // cout<<"sp_area="<<sp_area<<endl;

        // writeVTKoh("overhung.vtk",mv,mf,sp.oh_fn);

        // writeVTK("oh_pair.vtk",mv,mf,oh_pair);

        // ExportVTK(G,NoRemovableFaces,"noremovalface.vtk");
    }

    cout<<"sp_area_min="<<sp_area_min<<"direction_min="<<min_vector<<endl;
    cout<<"sp_area_max="<<sp_area_max<<"direction_max="<<max_vector<<endl;
    writeVTKVector("direction.vtk",{0,0,0},min_vector);
    writeVTK("oh_pair.vtk",mv,mf,max_overhung);

    // 回転角の計算
    double yaw_min=atan2(min_vector.y(),min_vector.x()); // z周り
    double pitch_min=atan2(sqrt(min_vector.x()*min_vector.x()+min_vector.y()*min_vector.y()),min_vector.z()); // y周り
    double roll_min=0; // x軸周り
    cout<<"min angle={"<<radToDeg(roll_min)<<","<<radToDeg(pitch_min)<<","<<radToDeg(yaw_min)<<"}\n";

    double yaw_max=atan2(max_vector.y(),max_vector.x()); // z周り
    double pitch_max=atan2(sqrt(max_vector.x()*max_vector.x()+max_vector.y()*max_vector.y()),max_vector.z()); // y周り
    double roll_max=0; // x軸周り
    cout<<"max angle={"<<radToDeg(roll_max)<<","<<radToDeg(pitch_max)<<","<<radToDeg(yaw_max)<<"}\n";

    ExportVTK(G,NoremovalFaces_min,"noremovalface_min.vtk");
    ExportVTK(G,NoremovalFaces_max,"noremocalface_max.vtk");
    
}