#include<vector>
#include<Eigen/Core>

class ro{
    public:

    std::vector<Eigen::Vector3d> vertex;
    std::vector<Eigen::Vector3i> face;
    std::string fn;

    void readPoint(std::vector<Eigen::Vector3d> &vertex,std::vector<Eigen::Vector3i> &face,std::string fn);
};