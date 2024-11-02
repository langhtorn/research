#include <Eigen/Dense>
#include <iostream>
#include <vector>

using namespace Eigen;
using namespace std;

// 点だけobjファイル出力
void obj_out(vector<Vector3d> vert,const char* file_name){
    FILE* f;
    f=fopen(file_name,"w");
    
    for(int i=0;i<vert.size();i++){
        // cout<<vert[i](0)<<","<<vert[i](1)<<endl;
        fprintf(f,"v %lf %lf %lf\n",vert[i](0),vert[i](1),vert[i](2));
    }
}

int main() {
    // 直方体の中心
    Vector3d center(1.0, 2.0, 3.0);

    // 高さを表す方向ベクトル
    Vector3d heightDir(1.0, 1.0, 1.0);  // 高さ方向

    // 高さに垂直な2つのベクトルを作る（幅と奥行きを表すベクトル）
    Vector3d widthDir(4.0, 0.0, 0.0);   // 幅方向
    Vector3d depthDir = heightDir.cross(widthDir).normalized() * 3.0; // 奥行き方向

    // 各ベクトルを0.5倍（中心からの距離）にして使う
    Vector3d halfHeight = 0.5 * heightDir;
    Vector3d halfWidth = 0.5 * widthDir;
    Vector3d halfDepth = 0.5 * depthDir;

    // 直方体の8つの頂点を計算
    vector<Vector3d> vertices;
    vertices.push_back(center + halfHeight + halfWidth + halfDepth); // + + +
    vertices.push_back(center + halfHeight + halfWidth - halfDepth); // + + -
    vertices.push_back(center + halfHeight - halfWidth + halfDepth); // + - +
    vertices.push_back(center + halfHeight - halfWidth - halfDepth); // + - -
    vertices.push_back(center - halfHeight + halfWidth + halfDepth); // - + +
    vertices.push_back(center - halfHeight + halfWidth - halfDepth); // - + -
    vertices.push_back(center - halfHeight - halfWidth + halfDepth); // - - +
    vertices.push_back(center - halfHeight - halfWidth - halfDepth); // - - -

    // 頂点出力
    for (int i = 0; i < vertices.size(); i++) {
        cout << "Vertex " << i + 1 << ": " << vertices[i].transpose() << endl;
    }

    obj_out(vertices,"example.obj");

    return 0;
}
