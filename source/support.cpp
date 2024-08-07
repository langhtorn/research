// サポート構造を構築する

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

using namespace std;
using namespace Eigen;

// 底面(x,z平面)でのグリッドセルを作る(入力：モデルの座標,グリッドの分割数)
vector<Vector3d> grid_cell(vector<Vector3d> V,int n,vector<Vector3d> &four){
    vector<Vector3d> gc;
    
    double x_min=100,x_max=-100,z_min=100,z_max=-100,y_min=100;
    for(int i=0;i<V.size();i++){
        if(x_min>V[i](0)){
            x_min=V[i](0);
        }
        if(z_min>V[i](2)){
           z_min=V[i](2); 
        }
        if(x_max<V[i](0)){
            x_max=V[i](0);
        }
        if(z_max<V[i](2)){
            z_max=V[i](2);
        }
        if(y_min>V[i](1)){
            y_min=V[i](1);
        }
    }
    four.push_back({x_min,y_min,z_min});
    four.push_back({x_max,y_min,z_min});
    four.push_back({x_min,y_min,z_max});
    four.push_back({x_max,y_min,z_max});

    Vector3d g;
    // cout<<"x_min="<<x_min<<" y_min="<<y_min<<endl;

    double x_length=x_max-x_min; //横の長さ
    double y_length=z_max-z_min; //縦の長さ
    for(double i=z_min;i<=z_max;i+=y_length/n){
        for(double j=x_min;j<=x_max;j+=x_length/n){
            // cout<<"("<<i<<" ,"<<j<<")"<<endl;
            g(0)=j;
            g(1)=y_min;
            g(2)=i;
            gc.push_back(g);
        }
    }

    // for(int i=0;i<V.size();i++){
    //     cout<<gc[i]<<endl;
    // }

    return gc;
}

// 点だけobjファイル出力
void obj_out(vector<Vector3d> vert,const char* file_name){
    FILE* f;
    f=fopen(file_name,"w");
    
    for(int i=0;i<vert.size();i++){
        // cout<<vert[i](0)<<","<<vert[i](1)<<endl;
        fprintf(f,"v %lf %lf %lf\n",vert[i](0),vert[i](1),vert[i](2));
    }
}

void obj_outf(vector<Vector3d> V,vector<Vector3i> F,const char* file_name){
    FILE* f;
    f=fopen(file_name,"w");
    for(int i=0;i<V.size();i++){
        fprintf(f,"v %lf %lf %lf\n",V[i](0),V[i](1),V[i](2));
    }
    for(int i=0;i<F.size();i++){
        fprintf(f,"f %d %d %d\n",F[i](0)+1,F[i](1)+1,F[i](2)+1);
    }
}

// オーバーハング点の判定(入力：グリッドセルの点，座標情報，オーバハング面の構成)
bool point_in_out(Vector3d P,vector<Vector3d> V,vector<Vector3i> oh_F,double h,vector<Vector3d> Vv){
    for(int i=0;i<oh_F.size();i++){

        Vector3d BA=V[oh_F[i](0)]-V[oh_F[i](1)];
        Vector3d BP=P-V[oh_F[i](1)];
        Vector3d BC=V[oh_F[i](2)]-V[oh_F[i](1)];
        Vector3d CP=P-V[oh_F[i](2)];
        Vector3d CA=V[oh_F[i](0)]-V[oh_F[i](2)];


        Vector3d cross_a=BP.cross(BA);
        // cout<<"BP*BA={"<<cross_a<<"}"<<endl;
        Vector3d cross_c=BC.cross(BP);
        // cout<<"BC*BP={"<<cross_c<<"}"<<endl;
        Vector3d cross_b=CA.cross(CP);
        // cout<<"CA*CP={"<<cross_b<<"}"<<endl;

        // cout<<"a="<<cross_a(2)<<" b="<<cross_b(2)<<" c="<<cross_c(2)<<endl;

        if((cross_a(1)>=0 && cross_b(1)>=0 && cross_c(1)>=0) || (cross_a(1)<=0 && cross_b(1)<=0 && cross_c(1)<=0)){ // 外積の向きが揃うとき内側
            // cout<<"内側\n";
            if(Vv[oh_F[i](0)](1)-P(1)>h*0.01) return true;
            // 面と地面の距離が閾値以下の場合，オーバーハング点にしない．サポートを立てる必要がないから．
        }

    }
    // cout<<"外側\n";
    return false;
}

// オーバハング面の面番号を返す
vector<int> oh_Fnum(double angle,Vector3d direction,vector<Vector3d> V,vector<Vector3i> F){

    vector<int> of;

    for(int i=0;i<F.size();i++){
        if(overh(angle,direction,V[F[i](0)],V[F[i](1)],V[F[i](2)])){ //オーバハング面だった場合
            of.push_back(i); //面番号
        }
    }

    return of;
}

// オーバーハング点の番号を返す(グリッドセルの何番の点か⇒↓)
vector<int> oh_Vnum(vector<Vector3d> gp,vector<Vector3d> V,vector<Vector3i> oh_F,double h,vector<Vector3d> Vv){
    vector<int> ohv;

    for(int i=0;i<gp.size();i++){
        if(point_in_out(gp[i],V,oh_F,h,Vv)) ohv.push_back(i); //オーバハング点だった場合
    }

    return ohv;
}

//  四隅の点の作成
vector<Vector3d>  f_corners(Vector3d P,int Nc,vector<Vector3d> four){
    double x_length=four[1](0)-four[0](0);
    double z_length=four[2](2)-four[0](2);
    double ox=(x_length/Nc)/2;
    double oz=(z_length/Nc)/2;

    vector<Vector3d> fc;
    fc.push_back({P(0)-ox,four[0](1),P(2)-oz});
    fc.push_back({P(0)+ox,four[0](1),P(2)-oz});
    fc.push_back({P(0)-ox,four[0](1),P(2)+oz});
    fc.push_back({P(0)+ox,four[0](1),P(2)+oz});

    return fc;
}

// 柱作成(柱の座標，柱の面,底面の頂点座標,柱の高さ,面番号始まりの点)
void pillar(vector<Vector3d> &pv,vector<Vector3i> &pf,vector<Vector3d> fc,double h,int ff){
    for(int i=0;i<4;i++){
        pv.push_back(fc[i]);
    }

    for(int i=0;i<4;i++){
        pv.push_back({fc[i](0),fc[i](1)+h,fc[i](2)});
    }

    pf.push_back({ff,ff+3,ff+2});
    pf.push_back({ff+3,ff+0,ff+1});
    pf.push_back({ff+6,ff+2,ff+7});
    pf.push_back({ff+7,ff+2,ff+3});
    pf.push_back({ff+7,ff+3,ff+5});
    pf.push_back({ff+5,ff+3,ff+1});
    pf.push_back({ff+5,ff+1,ff+4});
    pf.push_back({ff+4,ff+1,ff});
    pf.push_back({ff+4,ff,ff+6});
    pf.push_back({ff+6,ff,ff+2});
    pf.push_back({ff+4,ff+6,ff+5});
    pf.push_back({ff+5,ff+6,ff+7});

}

// サポート構築(Lattice) 入力：オーバーハング点，AABBの底面，セルの分割数，柱の高さ，求めたいサポートの点と面
void L_support(vector<Vector3d> ohp,vector<Vector3d> four,int Nc,double h,vector<Vector3d> &support_v,vector<Vector3i> &support_f,vector<Vector3d> V,vector<Vector3i> F){

    int f_num=0;
    for(int i=0;i<ohp.size();i++){
        vector<Vector3d> fc=f_corners(ohp[i],Nc,four); // 四隅の点を求める

        vector<Vector3d> pv;
        vector<Vector3i> pf;
        pillar(pv,pf,fc,h,f_num);
        // サポートの配列に一本のサポート柱の情報を入れる
        support_v.insert(support_v.end(),pv.begin(),pv.end());
        support_f.insert(support_f.end(),pf.begin(),pf.end());
        // 配列の中身削除
        pv.clear();
        pf.clear();
        // 全オーバーハング点で繰り返し柱作成
        f_num+=8; //一つの柱に12面あるから
    }

    MatrixXd spv;
    MatrixXi spf;
    vt vtom;
    vtom.vmat(support_v,support_f,spv,spf);
    cout<<"サポート変換\n";
    igl::writeOBJ("support_L.obj",spv,spf);
    MatrixXd v;
    MatrixXi f;
    vtom.vmat(V,F,v,f);
    cout<<"元モデル変換\n";
    MatrixXd boolsupport_v;
    MatrixXi boolsupport_f;
    igl::MeshBooleanType boolean_type=static_cast<igl::MeshBooleanType>(2);
    cout<<"booleanタイプ選び\n";
    igl::copyleft::cgal::mesh_boolean(spv,spf,v,f,boolean_type,boolsupport_v,boolsupport_f);
    igl::writeOBJ("boolsupport.obj",boolsupport_v,boolsupport_f);
}

// サポート構築(Tree)
vector<Vector3d> T_support(){

}

int main(int argc,char* argv[])
{   
    string fn;
    FILE *fp;
    vector<Vector3d> V;
    vector<Vector3i> F;
    ro reado;

    if(argc==1){
        cout<<"入力してください\n";
        exit(1);
    }else{
        fn=argv[1];
        fp=fopen(argv[1],"r");
    }
    reado.readPoint(V,F,fn);
    cout<<"読み込み完了\n";
    // for(int i=0;i<V.size();i++){
    //     cout<<V[i]<<endl;
    // }
    // for(int i=0;i<F.size();i++){
    //     cout<<F[i]<<endl;
    // }

    // x,y平面グリッドセルの作成
    int N_cell=30;
    vector<Vector3d> four;
    vector<Vector3d> gc=grid_cell(V,N_cell,four);
    cout<<"グリッドセルの作成完了\n";
    // obj_out(gc,"grid_cell.obj"); //obj出力

    double y_min=100;
    double y_max=-100;
    for(int i=0;i<V.size();i++){ //y座標最小値
        if(y_min>V[i](1)){
            y_min=V[i](1);
        }
        if(y_max<V[i](1)){
            y_max=V[i](1);
        }
    }

    // 面のオーバハングを調べる
    cout<<"面のオーバーハングを調べる"<<F.size()<<endl;
    double angle=0.6108652381980153; //閾値:cura35度(0.6108652381980153)，45度：0.7853981633974483
    Vector3d direction(0,-1,0); //造形方向ベクトル
    vector<int> num; //その面がオーバーハングかどうか(1 or 0)
    for(int i=0;i<F.size();i++){
        if(overh(angle,direction,V[F[i](0)],V[F[i](1)],V[F[i](2)])) num.push_back(1);
        else num.push_back(0);
    }

    // オーバハング点の判定
    vector<int> oF=oh_Fnum(angle,direction,V,F); // オーバーハング面の集合を作る
    // Vの射影
    vector<Vector3d> Vs;
    for(int i=0;i<V.size();i++){
        Vector3d vss={V[i](0),y_min,V[i](2)};
        Vs.push_back(vss);
    }
    vector<Vector3i> oh_F; //オーバハングだけの面
    for(int i=0;i<oF.size();i++){
        Vector3i ohf={F[oF[i]](0),F[oF[i]](1),F[oF[i]](2)};
        oh_F.push_back(ohf);
    }
    double h=y_max-y_min;
    vector<int> oh_point=oh_Vnum(gc,Vs,oh_F,h,V); //オーバーハング点の点番号
    
    // デバッグ
    vector<Vector3d> ohp;
    cout<<"gc="<<gc.size()<<" ohp="<<oh_point.size()<<endl;
    for(int i=0;i<oh_point.size();i++){
        Vector3d op=gc[oh_point[i]];
        ohp.push_back(op);
    }

    obj_out(ohp,"oh_point.obj");

    vector<Vector3d> sv;
    vector<Vector3i> sf;
    cout<<"h="<<h<<endl;
    L_support(ohp,four,N_cell,h,sv,sf,V,F);
    cout<<"sv="<<sv.size()<<" sf="<<sf.size()<<endl;
    obj_outf(sv,sf,"support_L.obj");

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

// 面のクラス分け
// オーバハングの面だけ残す
// 隣の面がオーバハングか探索する
// https://docs.google.com/spreadsheets/d/1N3K9jsz0g74o5LToadMPtwWfvS1vS7ZjCbrNMR23uYA/edit?usp=sharing