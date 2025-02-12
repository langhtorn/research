#include<iostream>
#include<vector>
#include<Eigen/Core>
#include<Eigen/Geometry>
#include<Eigen/Dense>
#include "support.hpp"

class judge{
    public:

    spt sp; //サポート情報
    std::vector<Eigen::Vector3d> samplingV; //サンプリング点
    std::vector<Eigen::Vector3i> samplingF; //サンプリングした際にできる面
    double r; //球の半径

    // 球面サンプリング
    //辺の中点を計算する関数
    Eigen::Vector3d findMidpoint(Eigen::Vector3d A,Eigen::Vector3d B);
    // 正二十面体のの初期解
    std::vector<Eigen::Vector3i> generate_icosahedron_faces();
    std::vector<Eigen::Vector3d> generate_icosahedron_vertices(const Eigen::Vector3d& center_point);
    // サブディビジョンする(サブディビジョン回数)
    void subdivide(int n);
    // モデルAABBツリーと三角形メッシュとの交差判定
    bool intersect_triangle(Eigen::MatrixXd TV,Eigen::MatrixXi TF);
    // 中心点から周囲の１点を求める(中心点と正規化された方向ベクトル)
    Eigen::MatrixXd four_point(Eigen::Vector3d centerpoint,Eigen::Vector3d direction);
    // 削除できるサポートを探す
    void delete_suppport();

    
};