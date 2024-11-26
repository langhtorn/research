// 元モデルGと疑似ツールの交差判定のプログラム
// 1本だけためしに

#include<iostream>
#include<vector>
#include<Eigen/Dense>
#include<Eigen/Core>
#include<igl/copyleft/cgal/intersect_other.h>
#include<igl/readOBJ.h>

using namespace std;
using namespace Eigen;

struct Model{
    MatrixXd V;
    MatrixXi F;
};

int main()
{
    Model G;
    igl::readOBJ("rect_cube1.obj",G.V,G.F);
    Model T;
    igl::readOBJ("rect_cube2.obj",T.V,T.F);

    MatrixXi IF;
    bool kari=true;

    bool intersects=igl::copyleft::cgal::intersect_other(G.V,G.F,T.V,T.F,kari,IF);

    if(intersects) cout<<"交差する\n";
    else cout<<"交差しない\n";

    return 0;
}