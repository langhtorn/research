#include<Eigen/Core>
#include<vector>

// class ob{
//     public:

//     Eigen::Vector3d o;
//     Eigen::Vector3d a;
//     Eigen::Vector3d b;
//     Eigen::Vector3d direction;
//     double angle;

    Eigen::Vector3d normal(Eigen::Vector3d o,Eigen::Vector3d a,Eigen::Vector3d b);

    bool overh(double angle,Eigen::Vector3d direction,Eigen::Vector3d o,Eigen::Vector3d a,Eigen::Vector3d b);


// };