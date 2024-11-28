#include<Eigen/Core>
#include<igl/signed_distance.h>
#include<igl/readOBJ.h>

using namespace std;
using namespace Eigen;

int main()
{
    MatrixXd P;
    MatrixXi kari;
    MatrixXd V;
    MatrixXi F;

    igl::readOBJ("icosahedron.obj",P,kari);
    igl::readOBJ("kitten_5408.obj",V,F);

    VectorXd S;
    VectorXi I;
    MatrixXd C;
    VectorXd N;

    igl::signed_distance(P,V,F,igl::SignedDistanceType::SIGNED_DISTANCE_TYPE_DEFAULT,S,I,C,N);
}