#include<iostream>
#include<vector>
#include<Eigen/Core>
#include<Eigen/Geometry>
#include<Eigen/Dense>
#include "overhang.hpp"
#include<cmath>
#include<igl/copyleft/cgal/mesh_boolean.h>

class spt{
    public:

    std::vector<Eigen::Vector3d> VG; //モデルの点情報
    std::vector<Eigen::Vector3i> FG; //モデルの面情報
    std::vector<Eigen::Vector3d> Vs; //底面に射影したモデルの点座標
    int N; //グリッドの分割数
    std::vector<Eigen::Vector3d> gc_V; //グリッドセルの点
    std::vector<int> oh_fn; //オーバーハング点の面番号
    std::vector<Eigen::Vector3i> oh_F; //オーバーハングだけの面
    std::vector<Eigen::Vector3d> ohp; //オーバーハング点の座標情報
    std::vector<int> ohvn; //オーバーハング点の番号
    std::vector<std::vector<Eigen::Vector3d>> rays_s; //モデルとrays_structureの交点座標
    std::vector<double> height; //各オーバーハング点の高さ
    std::vector<Eigen::Vector3d> four; //モデルのmin.max
    Eigen::MatrixXd VG_Mver; //iglで使えるようにモデルの点座標のMatrixバージョンも作っておく
    Eigen::MatrixXi FG_Mver; //モデルの面リストのMatrixバージョン
    Eigen::MatrixXd SV; //サポート頂点Matrixバージョン
    Eigen::MatrixXi SF; //サポート面リストMatrixバージョン
    Eigen::MatrixXd b_v;
    Eigen::MatrixXi b_f;

    // 最小包含立体
    void aabb();

    // 底面(x,z平面)でのグリッドセルを作る(入力：モデルの座標,グリッドの分割数)
    void grid_cell();

    // 点だけobjファイル出力
    void obj_out(std::vector<Eigen::Vector3d> vert,const char* file_name);

    void obj_outf(std::vector<Eigen::Vector3d> V,std::vector<Eigen::Vector3i> F,const char* file_name);

    // それぞれの柱の高さを設定
    void l_height();

    // オーバーハング点の判定(入力：グリッドセルの点，座標情報，オーバハング面の構成)
    bool point_in_out(Eigen::Vector3d P,double h);

    // オーバハング面の面番号を返す
    void oh_Fnum(double angle,Eigen::Vector3d direction);

    // オーバーハング点の番号を返す(グリッドセルの何番の点か⇒↓)
    std::vector<int> oh_Vnum(std::vector<Eigen::Vector3d> gp,double h);

    //  四隅の点の作成
    std::vector<Eigen::Vector3d>  f_corners(Eigen::Vector3d P);

    // 柱作成(柱の座標，柱の面,底面の頂点座標,柱の高さ,面番号始まりの点)
    void pillar(std::vector<Eigen::Vector3d> &pv,std::vector<Eigen::Vector3i> &pf,std::vector<Eigen::Vector3d> fc,double h,int ff);

    // サポート構築(Lattice) 入力：オーバーハング点，AABBの底面，セルの分割数，柱の高さ，求めたいサポートの点と面
    void L_support();

};