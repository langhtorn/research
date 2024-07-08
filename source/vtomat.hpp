#include<vector>
#include<Eigen/Core>

class vt{
    public:

    std::vector<Eigen::Vector3d> A;
    std::vector<Eigen::Vector3i> B;
    Eigen::Matrix3d AA;
    Eigen::Matrix3i BB;

    void vmat(std::vector<Eigen::Vector3d> A,std::vector<Eigen::Vector3i> B,Eigen::Matrix3d &AA,Eigen::Matrix3i &BB);
};