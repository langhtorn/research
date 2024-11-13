#include<iostream>
#include<vector>
#include<Eigen/Core>
#include<Eigen/Geometry>
#include<Eigen/Dense>
#include "support.hpp"

class judge{
    public:

    spt sp; //サポート情報

    // 球面サンプリング
    std::vector<Eigen::Vector3d> SphericalSampling(Eigen::Vector3d p0,Eigen::Vector3d direction,Eigen::Vector3d r,int samples);
    // モデルAABBツリーと三角形メッシュとの交差判定
    bool intersect_triangle(Eigen::MatrixXd TV,Eigen::MatrixXi TF);
    // 中心点から周囲の１点を求める(中心点と正規化された方向ベクトル)
    Eigen::MatrixXd four_point(Eigen::Vector3d centerpoint,Eigen::Vector3d direction);
    // 削除できるサポートを探す
    void delete_suppport();

    
};