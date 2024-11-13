// vector<VectorXd>をMatrixXdに変換する

#include<iostream>
#include<Eigen/Core>
#include<vector>
#include"vtomat.hpp"

using namespace std;
using namespace Eigen;

void vt::vmat(std::vector<Eigen::Vector3d> A,std::vector<Eigen::Vector3i> B,Eigen::MatrixXd &AA,Eigen::MatrixXi &BB)
{
    MatrixXd ad(A.size(),3);
    MatrixXi bd(B.size(),3);
    for(int i=0;i<A.size();i++){
        for(int j=0;j<3;j++){
            ad(i,j)=A[i](j);
        }
    }
    for(int i=0;i<B.size();i++){
        for(int j=0;j<3;j++){
            bd(i,j)=B[i](j);
        }
    }
    AA=ad;
    BB=bd;
}

// int main()
// {
//     vector<Vector3d> a;
//     a.push_back({1,2,3});
//     a.push_back({3,4,6});
//     a.push_back({4,5,6});

//     vector<Vector3i> d;
//     d.push_back({1,2,3});
//     d.push_back({3,4,6});
//     d.push_back({4,5,6});
//     Matrix3d aa;
//     Matrix3i dd;

//     vmat(a,d,aa,dd);
//     cout<<"mat="<<aa<<endl;
// }