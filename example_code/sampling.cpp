// 正二十面体のサブディビジョン&疑似ツール生やしてみる

#include<Eigen/Dense>
#include<vector>
#include<map>
#include<cmath>
#include<iostream>
#include<fstream>
#include<string>
#include<igl/signed_distance.h>
#include<igl/readOBJ.h>

using namespace Eigen;
using namespace std;

struct Model{
    vector<Vector3d> V;
    vector<Vector3i> F;
    double r;
    vector<Vector3d> newV;
    vector<Vector3i> newF;
    MatrixXd GV;
    MatrixXi GF;

    // 辺の中点を計算する&球面上に配置
    Vector3d findMindpoint(Vector3d A,Vector3d B){
        Vector3d midpoint=(A+B)/2.0;
        return midpoint.normalized()*r; //球面上に正規化してスケーリング
    }

    // サブディビジョンする 引数：サブディビジョン回数
    void subdivide(int n){
        for(int j=0;j<n;j++){
            newV.clear();
            newF.clear();
            for(int i=0;i<F.size();i++){
                // cout<<i<<endl;
                int v1=F[i](0);
                int v2=F[i](1);
                int v3=F[i](2);

                // 中点頂点計算
                Vector3d m12v=findMindpoint(V[v1],V[v2]);
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
        MatrixXd nv(V.size(),3);
        for(int i=0;i<V.size();i++){
            nv.row(i)=V[i];
        } 
        MatrixXd GV_3D(GV.rows(), 3); // 3次元座標行列として再定義
        for (int i = 0; i < GV.rows(); i++) {
            GV_3D(i, 0) = GV(i, 0);
            GV_3D(i, 1) = GV(i, 1);
            GV_3D(i, 2) = GV(i, 2);
        }
        MatrixXi GF_3D(GF.rows(), 3); // 3つの頂点を持つ面行列として再定義
        for (int i = 0; i < GF.rows(); i++) {
            GF_3D(i, 0) = GF(i, 0);
            GF_3D(i, 1) = GF(i, 1);
            GF_3D(i, 2) = GF(i, 2);
        }
        // 内外判定して削除する
        VectorXd S; //符号付距離
        VectorXi I; //最近接面のインデックス
        MatrixXd C; //最近接点
        VectorXd N; //法線方向
        igl::signed_distance(nv,GV_3D,GF_3D,igl::SignedDistanceType::SIGNED_DISTANCE_TYPE_DEFAULT,S,I,C,N);
        for(int i=0;i<V.size();i++){
            if(S[i]<0) newF.erase(newF.begin()+i);
        }   
    }

};

std::vector<Eigen::Vector3i> generate_icosahedron_faces() {
    return {
        { 0, 11, 5 }, { 0, 5, 1 }, { 0, 1, 7 }, { 0, 7, 10 }, { 0, 10, 11 },
        { 1, 5, 9 }, { 5, 11, 4 }, { 11, 10, 2 }, { 10, 7, 6 }, { 7, 1, 8 },
        { 3, 9, 4 }, { 3, 4, 2 }, { 3, 2, 6 }, { 3, 6, 8 }, { 3, 8, 9 },
        { 4, 9, 5 }, { 2, 4, 11 }, { 6, 2, 10 }, { 8, 6, 7 }, { 9, 8, 1 }
    };
}

std::vector<Eigen::Vector3d> generate_icosahedron_vertices(const Eigen::Vector3d& center_point, double r) {
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

struct tool{
    vector<Vector3d> tV;
    vector<Vector3i> tF;

    // 中心点から周囲の１点を求める(中心点と正規化された方向ベクトル)
    vector<Vector3d> four_point(Vector3d centerpoint,Vector3d direction){

        Vector3d v0=direction.unitOrthogonal();
        Vector3d v1=direction.cross(v0).normalized(); //長さを等しくする
        vector<Vector3d> p;

        // 周囲の4点求める
        double scale=0.05; //中心点からの距離 横幅
        p.push_back(centerpoint+scale*v0+scale*v1);
        p.push_back(centerpoint+scale*v0-scale*v1);
        p.push_back(centerpoint-scale*v0-scale*v1);
        p.push_back(centerpoint-scale*v0+scale*v1);
        
        return p;
    }

    // 疑似ツールの作成(作りたい方向ベクトル,始点)
    void make_toll(Vector3d direction,Vector3d centerpoint,double L){
        direction.normalize();
        tV=four_point(centerpoint,direction);

        centerpoint=centerpoint+L*direction;
        vector<Vector3d> tp2=four_point(centerpoint,direction);
        tV.insert(tV.end(),tp2.begin(),tp2.end());

        // 面の作成
        tF.push_back({0,1,2});
        tF.push_back({0,2,3});
        tF.push_back({4,5,6});
        tF.push_back({4,6,7});
        tF.push_back({2,1,5});
        tF.push_back({2,5,6});
        tF.push_back({3,2,6});
        tF.push_back({3,6,7});
        tF.push_back({0,3,7});
        tF.push_back({0,7,4});
        tF.push_back({1,0,4});
        tF.push_back({1,4,5});
    }
};

int main()
{
    Model G;
    G.r=1.0;
    Vector3d centerpoint={0,0,0};
    G.V=generate_icosahedron_vertices(centerpoint,G.r);
    G.F=generate_icosahedron_faces();
    G.newV=G.V;
    G.newF=G.F;

    igl::readOBJ("kitten_5508.obj",G.GV,G.GF);

    // FILE* f;
    // f=fopen("icosahedron.obj","w");
    // cout<<G.V.size()<<endl;
    // for(int i=0;i<G.V.size();i++){
    //     fprintf(f,"v %lf %lf %lf\n",G.V[i](0),G.V[i](1),G.V[i](2));
    // }
    // cout<<G.F.size()<<endl;
    // for(int i=0;i<G.F.size();i++){
    //     fprintf(f,"f %d %d %d\n",G.F[i](0)+1,G.F[i](1)+1,G.F[i](2)+1);
    // }
    // fclose(f);

    G.subdivide(2);

    FILE* f;
    f=fopen("icosahedron.obj","w");
    for(int i=0;i<G.newV.size();i++){
        fprintf(f,"v %lf %lf %lf\n",G.newV[i](0),G.newV[i](1),G.newV[i](2));
    }
    for(int i=0;i<G.newF.size();i++){
        fprintf(f,"f %d %d %d\n",G.newF[i](0)+1,G.newF[i](1)+1,G.newF[i](2)+1);
    }
    fclose(f);

    // 疑似ツールの作成
    vector<tool> TL;

    for(int i=0;i<G.newV.size();i++){
        tool onetool;
        onetool.make_toll(G.newV[i],centerpoint,G.r);
        TL.push_back(onetool);
    }

    // 疑似ツールの出力
    FILE* f2;
    f2=fopen("karitool.obj","w");
    for(int i=0;i<TL.size();i++){
        // cout<<i<<endl;
        // cout<<TL[i].tV.size()<<endl;
        for(int j=0;j<TL[i].tV.size();j++){
            fprintf(f2,"v %lf %lf %lf\n",TL[i].tV[j](0),TL[i].tV[j](1),TL[i].tV[j](2));
        }       
    }
    int cnt=0; //面を増やさなければならない
    for(int i=0;i<TL.size();i++){
        // cout<<TL[i].tV.size()<<endl;
        for(int j=0;j<TL[i].tF.size();j++){
            fprintf(f2,"f %d %d %d\n",TL[i].tF[j](0)+1+cnt,TL[i].tF[j](1)+1+cnt,TL[i].tF[j](2)+1+cnt);
        }
        cnt+=8;
    }
    fclose(f2);

    return 0;
}