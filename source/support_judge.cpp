//サポートの良否判定

#include<iostream>
#include<random>
#include<vector>
#include<Eigen/Core>
#include<Eigen/Geometry>
#include<Eigen/Dense>
#include<cmath>
#include<igl/copyleft/cgal/mesh_boolean.h>
#include "readobj.hpp"
#include "overhang.hpp"
#include "vtomat.hpp"
#include<igl/writeOBJ.h>
#include "support.hpp"
#include "support_judge.hpp"
#include<igl/copyleft/cgal/intersect_other.h>

using namespace std;
using namespace Eigen;

// 球面サンプリング
vector<Vector3d> judge::SphericalSampling(Vector3d p0,Vector3d direction,Vector3d r,int samples){
    vector<Vector3d> sample_p;
    double ts=M_PI/samples; //θの間隔
    double ps=M_PI/samples; //φの間隔

    for(int i=0;i<samples;i++){
        double theta=i*ts;
        for(int j=0;j<samples;j++){
            double phi=j*ps;
            Vector3d sp;
            sp<<r*sin(theta)*cos(phi),r*sin(theta)*sin(phi),r*cos(theta);

            // サンプリング点を方向ベクトルに移動
            Vector3d sd=direction*r+sp;
            sample_p.push_back(p0+sd);
        }
    }
    return sample_p;
}

// モデルAABBツリーと三角形メッシュとの交差判定
bool judge::intersect_triangle(MatrixXd TV,MatrixXi TF){

    MatrixXi IF; // 交差する三角形のペアが格納される

    bool intersects=igl::copyleft::cgal::intersect_other(sp.VG_Mver,sp.FG_Mver,TV,TF,IF);

    if(intersects) return true;
    else return false;

}

// 中心点から周囲の１点を求める(中心点と正規化された方向ベクトル)
MatrixXd judge::four_point(Vector3d centerpoint,Vector3d direction){

    MatrixXd fp(1,3);
    Vector3d v0=direction.unitOrthogonal();
    Vector3d v1=direction.cross(v0);
    vector<Vector3d> p;

    // 周囲の4点求める
    double scale=1; //中心点からの距離
    p.push_back(centerpoint+scale*v0+scale*v1);
    p.push_back(centerpoint+scale*v0-scale*v1);
    p.push_back(centerpoint-scale*v0+scale*v1);
    p.push_back(centerpoint-scale*v0-scale*v1);

    fp.row(0)=p[0].transpose();

    for(int i=1;i<4;i++){
        fp.conservativeResize(fp.rows()+1,fp.cols());
        fp.row(fp.rows()-1)=p[i].transpose();
    }
    
    return fp;
}

// 削除できるサポートを探す
void judge::delete_suppport(){
    for(int i=0;i<sp.ohvn.size();i++){ //オーバーハング点の番号，サポート番号

        for(int j=0;j<sp.rays_s[i].size();j++){ //始点を決める
            // 道具ボクセルの作成
            
            Vector3d direction(0,-1,0); // 初期方向ベクトルは造形方向下向きで
            Vector3d centerpoint=sp.rays_s[i][j]; //始点を平面の中心点とする

            direction.normalize(); //正規化して単位ベクトルにする
            
            MatrixXd toolpoints; //ツールの頂点座標
            MatrixXi toolface; //ツールの面座標

            toolpoints=four_point(centerpoint,direction); //周囲の4点を求める

            double L=1; //ツールの長さ
            centerpoint=centerpoint+L*direction; // ツールボクセルの反対側

            bool intersects=intersect_triangle(toolpoints,toolface); //交差判定 trueが交差してる
        }
    }
}

int main(int argc,char* argv[])
{   
    string fn;
    FILE *fp;
    vector<Vector3d> V;
    vector<Vector3i> F;
    ro reado;
    spt sp;

    if(argc==1){
        cout<<"入力してください\n";
        exit(1);
    }else{
        fn=argv[1];
        fp=fopen(argv[1],"r");
    }
    reado.readPoint(V,F,fn);
    sp.VG=V;
    sp.FG=F;
    cout<<"読み込み完了\n";
    // for(int i=0;i<V.size();i++){
    //     cout<<V[i]<<endl;
    // }
    // for(int i=0;i<F.size();i++){
    //     cout<<F[i]<<endl;
    // }

    // x,y平面グリッドセルの作成
    sp.N=30; //グリッドの分割数
    sp.aabb(); //最小包含
    sp.grid_cell();
    cout<<"グリッドセルの作成完了\n";
    // obj_out(gc,"grid_cell.obj"); //obj出力

    double y_max=-100;
    for(int i=0;i<V.size();i++){ //y座標最大
        if(y_max<V[i](1)){
            y_max=V[i](1);
        }
    }

    // 面のオーバハングを調べる
    cout<<"面のオーバーハングを調べる"<<sp.FG.size()<<endl;
    double angle=0.5235987755983; //閾値:cura35度(0.6108652381980153)，45度：0.7853981633974483,30度0.5235987755983
    Vector3d direction(0,-1,0); //造形方向ベクトル
    vector<int> num; //その面がオーバーハングかどうか(1 or 0)
    for(int i=0;i<F.size();i++){
        if(overh(angle,direction,V[F[i](0)],V[F[i](1)],V[F[i](2)])) num.push_back(1);
        else num.push_back(0);
    }

    // オーバハング点の判定
    sp.oh_Fnum(angle,direction); // オーバーハング面の集合を作る
    // Vの射影
    for(int i=0;i<V.size();i++){
        Vector3d vss={V[i](0),sp.four[0](1),V[i](2)};
        sp.Vs.push_back(vss);
    }

    for(int i=0;i<sp.oh_fn.size();i++){
        Vector3i ohf={F[sp.oh_fn[i]](0),F[sp.oh_fn[i]](1),F[sp.oh_fn[i]](2)};
        sp.oh_F.push_back(ohf);
    }
    double h=y_max-sp.four[0](1);
    vector<int> oh_point=sp.oh_Vnum(sp.gc_V,h); //オーバーハング点の点番号
    
    // デバッグ
    cout<<"gc="<<sp.gc_V.size()<<" ohp="<<oh_point.size()<<endl;
    for(int i=0;i<oh_point.size();i++){
        Vector3d op=sp.gc_V[oh_point[i]];
        sp.ohp.push_back(op);
    }
    cout<<"ohpsize="<<sp.ohp.size()<<endl;

    sp.obj_out(sp.ohp,"oh_point.obj");

    vector<Vector3d> sv; //サポート点情報
    vector<Vector3i> sf; //サポート面情報
    sp.l_height();
    cout<<"h="<<h<<endl;
    sp.L_support();
    cout<<"sv="<<sv.size()<<" sf="<<sf.size()<<endl;
    sp.obj_outf(sv,sf,"support_L.obj");

    // // vtk出力，結果の可視化
    // FILE* fpv;
    // fpv=fopen("overhang.vtk","w");
    // fprintf(fpv,"# vtk DataFile Version 2.0\n");
    // fprintf(fpv,"Title of my super cool VTK file\nASCII\nDATASET UNSTRUCTURED_GRID\n");
    // fprintf(fpv,"POINTS %d float\n",(int)V.size());
    // for(int i=0;i<V.size();i++){
    //     fprintf(fpv,"%lf %lf %lf\n",V[i](0),V[i](1),V[i](2));
    // }
    // fprintf(fpv,"CELLS %d %d\n",(int)F.size(),(int)F.size()*4);
    // for(int i=0;i<F.size();i++){
    //     fprintf(fpv,"3 %d %d %d\n",F[i](0),F[i](1),F[i](2));
    // }
    // fprintf(fpv,"CELL_TYPES %d\n",(int)F.size());
    // for(int i=0;i<F.size();i++){
    //     fprintf(fpv,"5\n");
    // }
    // fprintf(fpv,"CELL_DATA %d\n",(int)F.size());
    // fprintf(fpv,"SCALARS my_face_attr int 1\nLOOKUP_TABLE default\n");
    // for(int i=0;i<num.size();i++){
    //     fprintf(fpv,"%d\n",num[i]);
    // }
    // cout<<"numsize="<<num.size()<<endl;
}