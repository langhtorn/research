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
    std::vector<Eigen::Vector3d> SphericalSampling(Vector3d p0,Vector3d direction,Vector3d r,int samples);
    // 削除できるサポートを探す
    void delete_suppport();

    
};