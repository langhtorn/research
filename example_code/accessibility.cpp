// アクセス可能性の計算

#include<iostream>
#include<vector>
#include<map>
#include<Eigen/Core>
#include<Eigen/Geometry>
#include<Eigen/Dense>
#include<cmath>
#include<igl/copyleft/cgal/convex_hull.h>
#include<boost/geometry.hpp>
#include<boost/geometry/index/rtree.hpp>
#include<boost/geometry/geometries/box.hpp>
#include<utility>
#include<igl/readOBJ.h>
#include<igl/readOFF.h>
#include<igl/writeOBJ.h>

using namespace Eigen;
using namespace std;

namespace bg=boost::geometry;
namespace bgi=boost::geometry::index;
typedef bg::model::point<double,2,bg::cs::cartesian> Point2D; // 2次元ポイント型(θ,Φ)
typedef bg::model::box<Point2D> Box2D; // 2次元の矩形型(θ,Φの範囲)
typedef pair<Point2D,int> RTreeValue; // R-treeの定義(Point2Dと元の頂点インデックス)
typedef bgi::rtree<RTreeValue,bgi::quadratic<16>> RTree;

// 球面長方形
struct Rectangle{
    double theta_min,theta_max; // θの範囲
    double phi_min,phi_max; // Φの範囲
    int kinds; // 矩形の種類(1:北極，2:南極，3:交差，4:通常)
};

// 元モデル
struct Model{
    vector<Vector3d> MV;
    vector<Vector3i> MF;
};

// 投影された点と元の頂点のインデックスを管理
struct ProjectPoint{
    vector<Vector3d> point; //投影された点
    vector<int> originalIndex; // 元モデルでのインデックス
    int Sphere_Index; // 単位球のインデックス
};

struct SphericalPolygon{
    vector<Vector3d> vertices; // 凸包の頂点
    MatrixXi cvface; // 凸包の面
};

struct InaccessRegion{
    int faceindex; //対象となる面fのインデックス
    SphericalPolygon region; // 面fのアクセス不可能領域I
};

// アクセシビリティ行列
using AccessStatusMatrix=vector<vector<bool>>;

struct AC{

    // クラス変数

    vector<Matrix<double,3,3>> S; // 単位球を分割した球面三角形のリスト(三点の座標を持つ行列)
    AccessStatusMatrix AccessStatus; // 球面三角形ごとのアクセシビリティを格納する行列(行:球面三角形の数 列:非凸面の数)
    double r; // 境界球の半径
    Model G; // 元モデルの情報 
    vector<InaccessRegion> I; // 面fに対してアクセス不可能な領域
    vector<vector<Rectangle>> rectangle_I; // 囲い込み球面矩形
    RTree rtree; //range tree
    vector<ProjectPoint> projectedPoints; //投影後の点
    int iI=0;
    int numiI=0;


    // デバック用objで可視化
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

    // vector<Matrix>の場合可視化obj
    void visualizeSphericalTrianglesToObj(const vector<Matrix<double, 3, 3>>& S, const string& filename) {
        ofstream objFile(filename);
    
        if (!objFile.is_open()) {
            cerr << "Failed to open file for writing: " << filename << endl;
            return;
        }

        // 頂点データの書き込み (v x y z)
        vector<Vector3d> vertices; // 頂点リストを作成
            for (const auto& triangle : S) {
                for (int i = 0; i < 3; ++i) {
                    Vector3d vertex = triangle.row(i); // 各三角形の3つの頂点
                    vertices.push_back(vertex);
                }
            }

            // 頂点を出力 (v x y z)
            for (const auto& vertex : vertices) {
                objFile << "v " << vertex.x() << " " << vertex.y() << " " << vertex.z() << endl;
            }

            // 面データの書き込み (f v1 v2 v3)
            int vertexOffset = 1;  // OBJのインデックスは1から始まる
            int faceIndex = 0;
            for (const auto& triangle : S) {
                objFile << "f " 
                        << faceIndex * 3 + vertexOffset << " "   // 1番目の頂点
                        << faceIndex * 3 + 1 + vertexOffset << " " // 2番目の頂点
                        << faceIndex * 3 + 2 + vertexOffset << endl; // 3番目の頂点
                ++faceIndex;
            }

            objFile.close();
            cout << "OBJ file written to: " << filename << endl;
    }

    void exportAccessStatusF0ToVTK(const vector<Matrix<double, 3, 3>>& S, const Model& G, const AccessStatusMatrix& AccessStatus, const string& filename) {
    ofstream vtkFile(filename);

    if (!vtkFile) {
        cerr << "Failed to open file: " << filename << endl;
        return;
    }

    vtkFile << "# vtk DataFile Version 3.0\n";
    vtkFile << "Access Status Visualization for F0\n";
    vtkFile << "ASCII\n";
    vtkFile << "DATASET POLYDATA\n";

    // 1. 頂点リスト
    vector<Vector3d> vertices;
    for (const auto& triangle : S) {
        for (int i = 0; i < 3; ++i) {
            vertices.push_back(triangle.row(i));
        }
    }

    vtkFile << "POINTS " << vertices.size() << " float\n";
    for (const auto& vertex : vertices) {
        vtkFile << vertex.x() << " " << vertex.y() << " " << vertex.z() << "\n";
    }

    // 2. 三角形のインデックスリスト
    vtkFile << "POLYGONS " << S.size() << " " << S.size() * 4 << "\n";
    for (size_t i = 0; i < S.size(); ++i) {
        vtkFile << "3 " << i * 3 << " " << i * 3 + 1 << " " << i * 3 + 2 << "\n";
    }

    // 3. スカラー値（色分け用）
    vtkFile << "CELL_DATA " << S.size() << "\n";
    vtkFile << "SCALARS access_status_F0 int 1\n"; // 整数値 (0: 赤, 1: 緑)
    vtkFile << "LOOKUP_TABLE default\n";

    // 面 F[0] を基準に AccessStatus[0][j] を元に色付け
    for (size_t j = 0; j < S.size(); ++j) {
        if (AccessStatus[j][0]) {
            vtkFile << "1\n"; // 緑
        } else {
            vtkFile << "0\n"; // 赤
        }
    }

    vtkFile.close();
    cout << "VTK file exported for F0: " << filename << endl;
}

    // MatrixXiをvector<Vector3d>に変換
    vector<Vector3d> matrixtovector(const MatrixXd &matrix){
        vector<Vector3d> result;

        for(int i=0;i<matrix.rows();i++){
            result.push_back(matrix.row(i));
        }
        return result;
    }
    // obj vector<Vector3d> MatrixXiの場合
    // 頂点と面をOBJファイルに書き出す関数
void writeToOBJ(const vector<Vector3d>& vertices, const MatrixXi& faces, const string& filename) {
    // 出力ファイルを開く
    ofstream file(filename);
    if (!file.is_open()) {
        cerr << "Failed to open file: " << filename << endl;
        return;
    }

    // 頂点の書き込み
    for (const auto& vertex : vertices) {
        file << "v " << vertex.x() << " " << vertex.y() << " " << vertex.z() << endl;
    }

    // 面の書き込み
    for (int i = 0; i < faces.rows(); i++) {
        file << "f " << faces(i, 0) + 1 << " " << faces(i, 1) + 1 << " " << faces(i, 2) + 1 << endl;
        // OBJのインデックスは1から始まるので、インデックスを+1して書き出します
    }

    // ファイルを閉じる
    file.close();
    cout << "OBJ file written to: " << filename << endl;
}

    // 球面矩形をメッシュ化してOBJファイルに保存する関数
    void saveRectangleAsOBJ(const Rectangle& rect,const string& filename,int resolution=10){
        vector<Vector3d> vertices; //頂点リスト
        vector<Vector3i> faces; //面リスト

        // θとΦに分割して頂点を生成
        // cout<<"θとΦに分割して頂点生成\n";
        for(int i=0;i<=resolution;i++){
            double theta=rect.theta_min + (rect.theta_max - rect.theta_min) * i / resolution;
            for(int j=0;j<=resolution;j++){
                double phi = rect.phi_min + (rect.phi_max - rect.phi_min) * j / resolution;

                // 球面座標をデカルト座標に変換
                double x = sin(theta) * cos(phi);
                double y = sin(theta) * sin(phi);
                double z = cos(theta);

                vertices.emplace_back(x, y, z);
            }
        }

        // 面の生成
        // cout<<"面の生成\n";
        for (int i = 0; i < resolution; ++i) {
            for (int j = 0; j < resolution; ++j) {
                int idx1 = i * (resolution + 1) + j;
                int idx2 = idx1 + 1;
                int idx3 = (i + 1) * (resolution + 1) + j;
                int idx4 = idx3 + 1;

                // 2つの三角形に分割
                faces.emplace_back(idx1, idx2, idx3);
                faces.emplace_back(idx2, idx4, idx3);
            }
        }

        // OBJファイルに書き込む
        // cout<<"Obj化\n";
        ofstream obj_file(filename);
        if (!obj_file) {
            cerr << "Error: Could not open file " << filename << " for writing." << endl;
            return;
        }

        // 頂点の書き込み
        // cout<<"頂点\n";
        for (const auto& v : vertices) {
            obj_file << "v " << v.x() << " " << v.y() << " " << v.z() << endl;
        }

        // 面の書き込み
        // cout<<"面\n";
        for (const auto& f : faces) {
            obj_file << "f " << (f.x() + 1) << " " << (f.y() + 1) << " " << (f.z() + 1) << endl;
        }

        obj_file.close();
        cout << "OBJ file saved: " << filename << endl;

    }

    // ----------------------------------------------------------
    // 1.ジオデシック分割
    // 面の作成
    std::vector<Eigen::Vector3i> generate_icosahedron_faces() {
        return {
            { 0, 11, 5 }, { 0, 5, 1 }, { 0, 1, 7 }, { 0, 7, 10 }, { 0, 10, 11 },
            { 1, 5, 9 }, { 5, 11, 4 }, { 11, 10, 2 }, { 10, 7, 6 }, { 7, 1, 8 },
            { 3, 9, 4 }, { 3, 4, 2 }, { 3, 2, 6 }, { 3, 6, 8 }, { 3, 8, 9 },
            { 4, 9, 5 }, { 2, 4, 11 }, { 6, 2, 10 }, { 8, 6, 7 }, { 9, 8, 1 }
        };
    }
    // 頂点の作成
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
    // 辺の中点を計算する&球面上に配置&newVまで作成
    int findMindpoint(int v1,int v2,const vector<Vector3d>& V,vector<Vector3d>& newV,map<pair<int,int>,int>& midpoint_map){

        // 頂点インデックスをソート
        pair<int, int> edge = (v1 < v2) ? make_pair(v1, v2) : make_pair(v2, v1);

        // 重複チェック
        if(midpoint_map.count(edge)) return midpoint_map[edge]; // すでに存在していた場合，中点のインデックスを返す

        // 中点の計算
        Vector3d midpoint=((V[v1]+V[v2])/2.0);
        midpoint=midpoint.normalized();

        // 新しい頂点の登録
        int midpoint_index=newV.size();
        newV.push_back(midpoint);

        // マップに登録
        midpoint_map[edge]=midpoint_index;

        return midpoint_index;

    }
    // 球面三角形作成(n:分割数)
    void make_SphericalTriangles(int n){

        // 初期二十面体の生成
        Vector3d center_point={0,0,0};
        vector<Vector3i> F=generate_icosahedron_faces();
        vector<Vector3d> V=generate_icosahedron_vertices(center_point);

        // visualizeMeshToObj(V,F,"syokimenn.obj");

        // 分割のループ
        for(int j=0;j<n;j++){
            vector<Vector3i> newF;
            vector<Vector3d> newV=V;
            map<pair<int,int>,int> midpoint_map; // 中点の重複チェック

            for(const auto& face : F){
                // cout<<i<<endl;
                int v1=face(0);
                int v2=face(1);
                int v3=face(2);

                // 中点頂点計算
                int m12=findMindpoint(v1,v2,V,newV,midpoint_map);
                int m23=findMindpoint(v2,v3,V,newV,midpoint_map);
                int m31=findMindpoint(v3,v1,V,newV,midpoint_map);

                // 新しい三角形を生成
                newF.push_back({v1,m12,m31});
                newF.push_back({v2,m23,m12});
                newF.push_back({v3,m31,m23});
                newF.push_back({m12,m23,m31});

            }
            F=newF;
            V=newV;
        }
        visualizeMeshToObj(V,F,"syokimenn.obj");


        // sに三角形の頂点情報を格納
        for(const auto& face : F){
            Matrix<double,3,3> triangle;
            triangle.row(0)=V[face[0]];
            triangle.row(1)=V[face[1]];
            triangle.row(2)=V[face[2]];
            S.push_back(triangle);
        }
    }

    // 3.AccessStatus行列の初期値作成
    // 面Fに対して法線ベクトルを求める
    Vector3d calculateNormal(Vector3i f){
        Vector3d AB=G.MV[f(1)]-G.MV[f(0)];
        Vector3d AC=G.MV[f(2)]-G.MV[f(0)];

        Vector3d normal=AB.cross(AC);
        normal.normalize();

        return normal;
    }
    // 3.1 半球H(f)を計算する関数
    bool isPointInHemisphire(const Vector3d P,const Vector3d& normal){
        // 法線ベクトルと点の内積が0以上であればその点は半球に含まれる
        return P.dot(normal)>=0;
    }
    // 3.2 各面fに対して
    void setAccesStatus(){
        for(int i=0;i<G.MF.size();i++){
            Vector3d normal=calculateNormal(G.MF[i]); //面fの法線ベクトル

            // 面fの法線ベクトルを極とする半球H(f)に基づいて，各球面三角形sに対して判定
            for(int j=0;j<S.size();j++){ // 各球面三角形s
                Matrix<double,3,3> triangle=S[j]; // 三角形の頂点情報

                // 三角形の頂点情報を取得
                Vector3d v1=triangle.row(0);
                Vector3d v2=triangle.row(1);
                Vector3d v3=triangle.row(2);

                // すべての頂点が半球H(f)に含まれるかどうかチェック
                bool allPointInHemisphire=isPointInHemisphire(v1,normal) && isPointInHemisphire(v2,normal) && isPointInHemisphire(v3,normal);

                // 完全に含まれていたらTRUE
                if(allPointInHemisphire) AccessStatus[j][i]=true;
                else AccessStatus[j][i]=false;
            }
        }
    }

    // 初期化(Fn:非凸面の法線ベクトル)
    void INITIALIZE(){

        // visualizeMeshToObj(G.MV,G.MF,"motomesh.obj");
        // 1.境界球の表面を，球面三角形に分割する
        make_SphericalTriangles(3);

        visualizeSphericalTrianglesToObj(S,"sphericaltriangle.obj");

        // 2.AccessStatus行列を作成する
        AccessStatus.resize(S.size(),vector<bool>(G.MF.size(),false));

        // 3.AccessStatus行列の初期値を設定する
        setAccesStatus();

        // exportAccessStatusF0ToVTK(S, G, AccessStatus, "access_status_F0.vtk");

    }

    
    //---------------------------------------------------------- 
    // アクセス不可能領域Iを計算するアルゴリズム
    // 単位球の構築
    vector<Vector3d> constructUnitSpheres(){
        vector<Vector3d> unitSphereCenters;

        // 各頂点を単位球の中心とする
        for(int i=0;i<G.MV.size();i++){
            unitSphereCenters.push_back(G.MV[i].normalized());
        }

        return unitSphereCenters;
    }
    // 単位球への投影,投影後の点を返す(投影したい元の点，単位球の中心点)
    Vector3d ProjectOntoSphere(const Vector3d& point,const Vector3d& sphereCenter){
        // 点を正規化して単位球上に投影
        Vector3d direction=(point-sphereCenter).normalized();
        return direction;
    }
    // 面f'を指定された中心点の単位球に投影する(引数：面f単位球の中心点)
    void projectFaceOntoSpheres(vector<Vector3d> unitSpheres){
        // cout<<"単位球への投影\n";
        // cout<<"unitSphereSize="<<unitSpheres.size()<<endl;
        for(int i=0;i<unitSpheres.size();i++){

            vector<Vector3d> pv;
            vector<int> original_i;

            for(int j=0;j<G.MV.size();j++){

                
                // 面Fの頂点は投影しない
                if(i==j){
                    // cout<<"Skipping projection for face"<<j<<"on unit sphere"<<i<<endl;
                    continue; // なにもしない
                }else{
                    pv.push_back(ProjectOntoSphere(G.MV[j],unitSpheres[i]));
                    original_i.push_back(j);
                }
            }

            projectedPoints.push_back(ProjectPoint{pv,original_i,i});
        }
    }
    // 凸包を計算する関数(引数：面fの単位球に投影された点　)
    void computeConvexHull(const int facenum,const vector<Vector3d>& pp){
        // cout<<"凸包の計算\n";
        // 投影された点をMatrixXdに変換
        MatrixXd points(pp.size(),3);
        for(size_t i=0;i<pp.size();i++){
            points.row(i)=pp[i];
        }
        // cout<<"Matrixへの変換\n";

        // 凸包の計算結果を格納する行列
        MatrixXi faces; // 凸包を構成する三角形のインデックス

        igl::copyleft::cgal::convex_hull(points,faces);

        // if(facenum==0) igl::writeOBJ("cvhull1.obj",points,faces);

        vector<Vector3d> cvv=matrixtovector(points);

        // f が f' によってアクセス不可能な方向の集合 I は、CV の内部領域として定義する
        I.push_back(InaccessRegion{facenum,SphericalPolygon{cvv,faces}});

    }
    // 交差点が三角形の内部にあるか確認
    bool isPointInsideTriangle(const Vector3d& P, const Vector3d& A, const Vector3d& B, const Vector3d& C) {
        // ベクトルを計算
        Vector3d v0 = C - A;
        Vector3d v1 = B - A;
        Vector3d v2 = P - A;

        // バリデーション座標で三角形の内部判定
        double dot00 = v0.dot(v0);
        double dot01 = v0.dot(v1);
        double dot02 = v0.dot(v2);
        double dot11 = v1.dot(v1);
        double dot12 = v1.dot(v2);

        double invDenom = 1.0 / (dot00 * dot11 - dot01 * dot01);
        double u = (dot11 * dot02 - dot01 * dot12) * invDenom;
        double v = (dot00 * dot12 - dot01 * dot02) * invDenom;

        // u, v, u + v がすべて [0, 1] の範囲内なら三角形内部
        return (u >= -1e-6) && (v >= -1e-6) && (u + v <= 1.0 + 1e-6);

    }

    // z軸と交差するか確認
    int dosePlaneIntersectZAxis(const vector<Vector3d>& CV,const MatrixXi& faces){
        for(int i=0;i<faces.rows();i++){
            Vector3d A=CV[faces(i,0)];
            Vector3d B=CV[faces(i,1)];
            Vector3d C=CV[faces(i,2)];

            Vector3d AB=B-A;
            Vector3d AC=C-A;
            Vector3d normal=AB.cross(AC);

            // z軸方向の法線ベクトルが0でないことを確認
            if(abs(normal.z())<1e-6){
                continue; // z軸と平行なので交差しない
            }

            // 平面の方程式
            double d=normal.x()*A.x()+normal.y()*A.y()+normal.z()*A.z();

            // z軸との交差点のz座標を計算
            double intersection_z=d/normal.z();

            // 交差点
            Vector3d intersection_point(0,0,intersection_z);

            if(isPointInsideTriangle(intersection_point,A,B,C)){
                if(intersection_z>1e-6){
                    return 1; // 北極
                }else if(intersection_z<1e-6){
                    return 2; //南極
                }
            }
        }
        return 3; // 交差しない
    }
    // 内部領域に含まれるか判定する関数
    bool isPointInsideConvexHull(const Vector3d& point,const vector<Vector3d>& CV,MatrixXi faces){
        // cout<<"col"<<faces.cols()<<"rows"<<faces.rows()<<endl;
        double epsilon=1e-8;
        for(int i=0;i<faces.rows();i++){
            Vector3d v1=CV[faces(i,0)];
            Vector3d v2=CV[faces(i,1)];
            Vector3d v3=CV[faces(i,2)];

            // 面の法線を取得
            Vector3d facenormal=(v2-v1).cross(v3-v1).normalized();

            // 点と法線の内積で内外を判定
            // if(centroid.dot(point)>0){
            //     return true; // 外側にある
            // }
        }
        // 全ての面に対して内側であれば内部
        return false;
    }

    // INACCESSIBLEのmain
    void INACCESSIBLE(){

        // 1.ファセットfの3つの頂点に単位球の構築
        vector<Vector3d> unitSpheres=constructUnitSpheres();

        // 2.元モデルの各面f'を単位球に投影する
        projectFaceOntoSpheres(unitSpheres);
        cout<<"単位球への投影完了\n";

        vector<Vector3i> fa={{0,0,0}};
        // visualizeMeshToObj(projectedPoints[0].point,fa,"projectpoint.obj");

        // 各モデルの各面fに対して
        for(int i=0;i<G.MF.size();i++){
            // cout<<"面F "<<i<<endl;

            // 3.求めた頂点群に凸包を求める
            // 面Fの3つの投影点を結合する
            vector<Vector3d> pp3;
            pp3.insert(pp3.end(),projectedPoints[G.MF[i](0)].point.begin(),projectedPoints[G.MF[i](0)].point.end());
            pp3.insert(pp3.end(),projectedPoints[G.MF[i](1)].point.begin(),projectedPoints[G.MF[i](1)].point.end());
            pp3.insert(pp3.end(),projectedPoints[G.MF[i](2)].point.begin(),projectedPoints[G.MF[i](2)].point.end());
            if(i==29) visualizeMeshToObj(pp3,fa,"pp3.obj");
            computeConvexHull(i,pp3); // 凸包の計算とIの定義
            // cout<<"凸包面求めた\n";
                
        }
        cout<<"凸包出力\n";

        // 凸包の頂点と面のチェック
        writeToOBJ(I[29].region.vertices,I[29].region.cvface,"cvhull.obj");

    }


    // ----------------------------------------------------------
    // アクセス不可能領域Iに対して，囲い込み球面矩形R0を計算する関数
    vector<Rectangle> calculateEnclosingRectangle(const InaccessRegion& region){

        vector<Rectangle> rectangles;

        const vector<Vector3d>& vertices=region.region.vertices; //凸包の頂点
        // if(iI==100) writeToOBJ(vertices,region.region.cvface,"rectanglecv.obj");

        // Iが空の場合は処理しない
        if(vertices.empty()){
            // cerr<<"Error: InaccessRegion vector I is empty."<<endl;
            return rectangles;
        }

        // θ，Φの初期化
        double theta_min=numeric_limits<double>::infinity();
        double theta_max=-numeric_limits<double>::infinity();
        double phi_min=numeric_limits<double>::infinity();
        double phi_max=-numeric_limits<double>::infinity();

        // 各頂点の球面座標を計算&北極点と南極点を調べる
        for(const auto& vertex : vertices){

            // 球面座標に変換
            double sphere_r=vertex.norm(); // 半径
            double theta=acos(vertex.z()/sphere_r); // 緯度θ
            double phi=atan2(vertex.y(),vertex.x()); // 緯度Φ

            // Φを[0,2π]の範囲に変換
            // if(phi<0) phi+=2*M_PI;

            // θ，Φの範囲を更新
            theta_min=min(theta_min,theta);
            theta_max=max(theta_max,theta);
            phi_min=min(phi_min,phi);
            phi_max=max(phi_max,phi);

        }

        // 1.北極[0,0]が領域Iに含まれるか判定
        bool north_included=false;
        Vector3d north_pole(0,0,1);
        if(dosePlaneIntersectZAxis(vertices,region.region.cvface)==1){
            rectangles.push_back(Rectangle{0,theta_max,0,2*M_PI,1});
            north_included=true;
        }

        // 2.南極[π,0]が領域Iに含まれるか判定
        bool south_included=false;
        Vector3d south_pole(0,0,-1);
        if(dosePlaneIntersectZAxis(vertices,region.region.cvface)==2){
            rectangles.push_back(Rectangle{theta_min,M_PI,0,2*M_PI,2});
            south_included=true;
        }

        // 3. Φ=0の弧が領域Iに交差するか判定 極を通らなった場合に確認
        bool phi_split=false;
        if(phi_min<0 && phi_max>0 && north_included==false && south_included==false){
            // 矩形を二つに分割
            rectangles.push_back(Rectangle{theta_min,theta_max,0,phi_max,3});
            rectangles.push_back(Rectangle{theta_min,theta_max,phi_min+2*M_PI,2*M_PI,3});
            phi_split=true;
        }

        // 4.通常ケース　単一の矩形を生成
        if(!north_included && !south_included && !phi_split){
            if(phi_min<0) phi_min+=2*M_PI;
            if(phi_max<0) phi_max+=2*M_PI;
            rectangles.push_back(Rectangle{theta_min,theta_max,phi_min,phi_max,4});
        }

        // if(iI==100){
        //     cout<<"north:"<<north_included<<" south:"<<south_included<<" phi:"<<phi_split<<endl;
        //     cout<<"不可能領域のインデックス:"<<region.faceindex<<endl;
        // }
        iI++;
        return rectangles;
        
    } 
    void CALCULATEENCLOSINGRECTANGLE(){
        // writeToOBJ(I[100].region.vertices,I[100].region.cvface,"rectanglecv.obj");
        for(const auto& region : I){
            rectangle_I.push_back(calculateEnclosingRectangle(region));
        }
        // cout<<"rectangle_size="<<rectangle_I.size()<<endl;
        cout<<" 囲い込み球面矩形R0の確認\n";

        saveRectangleAsOBJ(rectangle_I[29][0],"enclosing_rectangle.obj");
    }


    // ----------------------------------------------------------
    // 占有テストを行う関数
    // パッチRの作成(元の矩形，拡張量)
    Rectangle generateCandidatePatch(const Rectangle& R0,double delta_theta,double delta_phi){
        Rectangle expanded;

        // θの範囲を拡張
        expanded.theta_min=R0.theta_min-delta_theta;
        expanded.theta_max=R0.theta_max+delta_theta;

        // Φの範囲を拡張
        if(R0.kinds==3 || R0.kinds==4){
            expanded.phi_min=R0.phi_min-delta_phi;
            expanded.phi_max=R0.phi_max+delta_phi;
        }else{
            expanded.phi_min=0;
            expanded.phi_max=2*M_PI;
        }

        // Φを[0,2π]の範囲に調整

        // if(expanded.phi_min<0){
        //     expanded.phi_min+=2*M_PI;
        //     karioki.phi_max=expanded.phi_min;
        // }
        // if(expanded.phi_max>2*M_PI){
        //     expanded.phi_max-=2*M_PI;
        //     karioki.phi_min=expanded.phi_max;
        // }
        // expanded.phi_min=karioki.phi_min;
        // expanded.phi_max=karioki.phi_max;

        return expanded;
    }
    // R-treeの構築
    void buildRTree(){
        int vnum=0;
        for(size_t i=0;i<S.size();i++){
            const auto& vertex=S[i];

            // 各三角形の3つの頂点を処理
            for(int j=0;j<3;j++){
                Vector3d vertex_i=vertex.row(j);
                // 球面座標に変換
                double sphere_r=vertex_i.norm(); // 半径
                double theta=acos(vertex_i.z()/sphere_r); // 緯度θ
                double phi=atan2(vertex_i.y(),vertex_i.x()); // 緯度Φ

                // Φを[0,2π]の範囲に変換
                if(phi<0) phi+=2*M_PI;

                // R-treeに(Point2D,元のインデックス)を挿入
                rtree.insert(make_pair(Point2D(theta,phi),vnum));
                vnum++;
            }

        }
    }
    // 範囲検索の実装
    vector<int> queryVerticesInRectangle(const Rectangle& R){
        vector<int> result_indices;


        // Rの範囲[θmin,θmax]x[Φmin,Φmax]をBox2Dに変換
        Box2D queryBox(Point2D(R.theta_min,R.phi_min),Point2D(R.theta_max,R.phi_max));

        // R-treeで範囲検索
        vector<RTreeValue> queryResults;
        rtree.query(bgi::within(queryBox),back_inserter(queryResults));

        // 検索結果からインデックスを抽出
        for(const auto& result:queryResults){
            result_indices.push_back(result.second);
        }

        return result_indices;
    }
    // 経度が0をまたぐ場合への対応
    vector<int> queryVerticesWithWrap(Rectangle& R){
        vector<int> result_indices;

        // 微小な誤差を許容する範囲を設定(これをしないとΦ=0の弧の点が含まれなくなる)
        const double epsilon = 1e-10;  // 誤差の許容範囲（適宜調整

        // Φの範囲がずれてたら[0,2π]に調整
        if(R.phi_min<0) R.phi_min+=2*M_PI;
        if(R.phi_max>2*M_PI) R.phi_max-=2*M_PI;

        // 極を含んだ場合は通常ケースで処理する
        if(R.kinds==1 || R.kinds==2){
            if(numiI==29){
                cout<<"通常\n";
                cout<<"R="<<R.phi_min<<","<<R.phi_max<<endl;
                saveRectangleAsOBJ(R,"R.obj");
            }
            // 通常ケースΦmin<=Φmax
            auto indices=queryVerticesInRectangle(R);
            result_indices.insert(result_indices.end(),indices.begin(),indices.end());
        }else{
            if(R.phi_min <= R.phi_max){
                // 通常ケースΦmin<=Φmax
                auto indices=queryVerticesInRectangle(R);
                result_indices.insert(result_indices.end(),indices.begin(),indices.end());
                
            }else{
                // Φmin>Φmaxの場合，二つの範囲に分割
                // cout<<"分割\n";
                Rectangle R1={R.theta_min,R.theta_max,0-epsilon,R.phi_max}; // θmin~2π
                Rectangle R2={R.theta_min,R.theta_max,R.phi_min,2*M_PI}; //0~Φmax

                auto indices1=queryVerticesInRectangle(R1);
                auto indices2=queryVerticesInRectangle(R2);
                if(numiI==29){
                    cout<<"分割\n";
                    saveRectangleAsOBJ(R1,"R1.obj");
                    cout<<"R1="<<R1.phi_min<<","<<R1.phi_max<<endl;
                    saveRectangleAsOBJ(R2,"R2.obj");
                    cout<<"R2="<<R2.phi_min<<","<<R2.phi_max<<endl;
                }

                // 結果を結合
                result_indices.insert(result_indices.end(),indices1.begin(),indices1.end());
                result_indices.insert(result_indices.end(),indices2.begin(),indices2.end());
            }
        }
        

        return result_indices;
    }
    // 球面三角形の面積を計算する関数
    double sphericalTriangleArea(const Vector3d& a,const Vector3d& b,const Vector3d& c){
        // 球面三角形の余弦定理に基づいて面積を計算
        double angleA=acos(b.dot(c)/(b.norm()*c.norm()));
        double angleB=acos(c.dot(a)/(c.norm()*a.norm()));
        double angleC=acos(a.dot(b)/(a.norm()*b.norm()));

        return angleA+angleB+angleC-M_PI;
    }
    // 球面三角形内に点が含まれるかの判定
    bool isPointSphericalTriangle(double theta,double phi,double theta1,double phi1,double theta2,double phi2,double theta3,double phi3){
        // 球面座標を直交座標系に変換
        auto sphericalToCartesian=[](double theta,double phi){
            return Vector3d(sin(theta)*cos(phi),sin(theta)*sin(phi),cos(theta));
        };

        Vector3d p=sphericalToCartesian(theta,phi);
        Vector3d a=sphericalToCartesian(theta1,phi1);
        Vector3d b=sphericalToCartesian(theta2,phi2);
        Vector3d c=sphericalToCartesian(theta3,phi3);

        // 三角形の面積の合計を計算
        double totalArea=sphericalTriangleArea(a,b,c);

        // 点を含む三角形の面積を計算
        double area1=sphericalTriangleArea(p,b,c);
        double area2=sphericalTriangleArea(a,p,c);
        double area3=sphericalTriangleArea(a,b,p);

        // 合計面積が元の三角形の面積と一致するか判定
        return abs(totalArea-(area1+area2+area3))<1e-6;
    }
    // 点が特定の領域面に含まれているかを判定する
    bool isWithinFace(double theta,double phi,const SphericalPolygon& region,const Vector3i& face){
        const auto& v1=region.vertices[face[0]];
        const auto& v2=region.vertices[face[1]];
        const auto& v3=region.vertices[face[2]];

        // 頂点の球面座標
        double r1=v1.norm();
        double theta1=acos(v1.z()/r1);
        double phi1=atan2(v1.y(),v1.x());
        if(phi1<0) phi1+=2*M_PI;

        double r2=v2.norm();
        double theta2=acos(v2.z()/r2);
        double phi2=atan2(v2.y(),v2.x());
        if(phi2<0) phi2+=2*M_PI;

        double r3=v3.norm();
        double theta3=acos(v3.z()/r3);
        double phi3=atan2(v3.y(),v3.x());
        if(phi3<0) phi3+=2*M_PI;

        // 球面三角形内に点が含まれているか確認
        return isPointSphericalTriangle(theta,phi,theta1,phi1,theta2,phi2,theta3,phi3);
    }
    // 点が領域に含まれるか判定
    bool isPointInRegion(const Vector3d& point,const SphericalPolygon& region){
        // 点の球面座標を計算
        double rr=point.norm();
        double theta=acos(point.z()/r); // 緯度θ
        double phi=atan2(point.y(),point.x()); // 経度Φ
        if(phi<0) phi+=2*M_PI; // Φを[0,2π]に正規化

        // θ，Φが領域の範囲に含まれるか判定
        for(int i=0;i<region.cvface.rows();i++){
            if(isWithinFace(theta,phi,region,region.cvface.row(i))){
                return true; // 含まれている
            }
        }
        return false; // 含まれていない
    }
    // アクセス可能性行列の更新と占有テスト(候補点VR,単位球の番号，面Fiのアクセス不可能領域)
    void updateAccessibilityMatrix(const vector<int>& VR,const int Sphere_i,const SphericalPolygon& inaccessibleRegion){
        // 点が領域内になるか判定する
        for(int i=0;i<VR.size();i++){
            // 三角形の頂点を取得

            int vertex_index=VR[i]; // 頂点インデックス

            // この頂点がどの三角形に属しているか探す
            for(int f_idx=0;f_idx<G.MF.size();f_idx++){
                const Vector3i& face=G.MF[f_idx]; // 三角形の頂点インデックス

                // 三角形の頂点インデックスの中にVR[i]が含まれていれば，この三角形に関連する
                if(face[0]==vertex_index || face[1]==vertex_index || face[2]==vertex_index){
                    const Vector3d& v0=projectedPoints[Sphere_i].point[face[0]];
                    const Vector3d& v1=projectedPoints[Sphere_i].point[face[1]];
                    const Vector3d& v2=projectedPoints[Sphere_i].point[face[2]];

                    // 頂点が全てアクセス不可能領域に含まれているか判定
                    bool isInaccessible=isPointInRegion(v0,inaccessibleRegion) && isPointInRegion(v1,inaccessibleRegion) && isPointInRegion(v2,inaccessibleRegion);

                    // 判定結果をAccessStatusに反映
                    if(isInaccessible){
                        // VR[i]が関連する三角形がアクセス不可能であれば，対応する行列を更新
                        AccessStatus[Sphere_i][f_idx]=true;
                    }
                    break;
                }
            }

        }
    }
    void OCCUPANCY(){

        // 1.アクセス不可能領域Iに対する囲い込み球面矩形R0を計算する
        // writeToOBJ(I[100].region.vertices,I[100].region.cvface,"occ.obj");
        CALCULATEENCLOSINGRECTANGLE();

        // 範囲木の構築
        buildRTree();

        // 解像度の計算
        double L=2*sin(M_PI/5);
        double n=sqrt(S.size()/20);
        double dL=L/n;

        
        for(int i=0;i<rectangle_I.size();i++){

            // 面Fiに対応するすべての囲い込み球面矩形を処理
            for(const Rectangle& R: rectangle_I[i]){

                // 2. 球面矩形R0を拡張して，候補パッチRを生成する
                Rectangle extendedRectangle=generateCandidatePatch(R,dL,dL);
                // cout<<"-------------------\n";
                if(numiI==29) saveRectangleAsOBJ(extendedRectangle,"expanded_Rectangle.obj");
                if(i%10==0) cout<<"候補パッチR"<<i<<" の作成\n";                


                // 3.range-treeを使用して，候補パッチRに含まれる球面頂点VRを取得する
                vector<int> VR=queryVerticesWithWrap(extendedRectangle);
                if(i%10==0) cout<<"VRsize="<<VR.size()<<endl;
                if(numiI==29){
                    // 頂点座標の取得
                    vector<Vector3d> VR_vert;
                    
                    for(const auto& index : VR){
                        // cout<<"VR:"<<index<<endl;
                        int index_i=index/3;
                        int index_j=index%3;

                        Vector3d vertex=S[index_i].row(index_j);
                        VR_vert.push_back(vertex);

                    }

                    vector<Vector3i> fa;
                    visualizeMeshToObj(VR_vert,fa,"VR_point.obj");
                }

                // 4.球面上の候補頂点および候補三角形に対して，占有テストを行う
                // 各面Fiのアクセス不可能領域に基づく領域を取得
                const SphericalPolygon& inaccessRegion=I[i].region;
                // 5.アクセス可能性行列の更新
                updateAccessibilityMatrix(VR,i,inaccessRegion);
                
            }
            numiI++;
        }

    }

    // アルゴリズム全般はここで
    void computeAccessibility(){

        // 初期化(球面三角形S,AccessStatus行列に関して)
        INITIALIZE();
        cout<<"初期化完了\n";

        // アクセス不可能領域を計算する関数
        INACCESSIBLE();
        cout<<"アクセス不可能領域計算完了\n";

        // 占有テストを行う関数
        OCCUPANCY();
        cout<<"占有テスト\n";
    }
};



int main()
{
    MatrixXd V;
    MatrixXi F;

    // ファイルの読み込み
    string filename="kitten_500.obj";
    if(!igl::readOBJ(filename,V,F)){
        cerr<<"Error\n"<<endl;
        return -1;
    }

    // Eigenからstd::vectorに変換
    vector<Vector3d> vertices; //元モデル頂点
    vector<Vector3i> faces; //元モデル面

    // 頂点の変換
    for(int i=0;i<V.rows();i++){
        vertices.emplace_back(V(i,0),V(i,1),V(i,2));
    }

    // 面の変換
    for(int i=0;i<F.rows();i++){
        faces.emplace_back(F(i,0),F(i,1),F(i,2));
    }
    cout<<"読み込み完了\n";

    Model originalModel;
    originalModel.MV=vertices;
    originalModel.MF=faces;

    // 球面分割情報の初期化
    vector<Matrix<double,3,3>> sphericalTriangles;
    
    // 投影点の初期化
    vector<ProjectPoint> projectedPoints;

    // AC構造体の初期化
    AC ac;
    ac.r=1.0; // 単位球の半径
    ac.G=originalModel;
    cout<<"初期化\n";

    // ac内のメインのアルゴリズム呼び出し
    ac.computeAccessibility();

    return 0;
}