// ミンコフスキー和をlibiglで実装

#include<Eigen/Core>
#include<Eigen/Geometry>
#include<igl/copyleft/cgal/minkowski_sum.h>
#include<igl/readOBJ.h>
#include<igl/writeOBJ.h>
#include<string>

using namespace std;
using namespace Eigen;

int main(int argc,char* argv[])
{
    string fn;
    FILE *fp;
    MatrixXd VA,VB;
    MatrixXi FA,FB;

    if(argc==1){
        cout<<"入力してください\n";
        exit(1);
    }else{
        igl::readOBJ(argv[1],VA,FA);
        igl::readOBJ(argv[2],VB,FB);
    }

    cout<<"読み込み完了\n";

    MatrixXd V;
    MatrixXi F;
    VectorXi J;
    igl::copyleft::cgal::minkowski_sum(VA,FA,VB,FB,true,V,F,J);

    cout<<"minkowski完了\n";

    igl::writeOBJ("minko.obj",V,F);

    return 0;
    
}