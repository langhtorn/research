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
#include "support.hpp"

using namespace std;
using namespace Eigen;

// 最小包含立体
void spt::aabb(){
    double x_min=100,x_max=-100,z_min=100,z_max=-100,y_min=100;
    for(int i=0;i<VG.size();i++){
        if(x_min>VG[i](0)){
            x_min=VG[i](0);
        }
        if(z_min>VG[i](2)){
           z_min=VG[i](2); 
        }
        if(x_max<VG[i](0)){
            x_max=VG[i](0);
        }
        if(z_max<VG[i](2)){
            z_max=VG[i](2);
        }
        if(y_min>VG[i](1)){
            y_min=VG[i](1);
        }
    }
    four.push_back({x_min,y_min,z_min});
    four.push_back({x_max,y_min,z_min});
    four.push_back({x_min,y_min,z_max});
    four.push_back({x_max,y_min,z_max});
}

// 底面(x,z平面)でのグリッドセルを作る(入力：モデルの座標,グリッドの分割数)
void spt::grid_cell(){

    Vector3d g;
    // cout<<"x_min="<<x_min<<" y_min="<<y_min<<endl;

    double x_length=four[3](0)-four[0](0); //横の長さ
    double y_length=four[3](2)-four[0](2); //縦の長さ
    for(double i=four[0](2);i<=four[3](2);i+=y_length/N){
        for(double j=four[0](0);j<=four[3](0);j+=x_length/N){
            // cout<<"("<<i<<" ,"<<j<<")"<<endl;
            g(0)=j;
            g(1)=four[0](1);
            g(2)=i;
            gc_V.push_back(g);
        }
    }
}

// 点だけobjファイル出力
void spt::obj_out(vector<Vector3d> vert,const char* file_name){
    FILE* f;
    f=fopen(file_name,"w");
    
    for(int i=0;i<vert.size();i++){
        // cout<<vert[i](0)<<","<<vert[i](1)<<endl;
        fprintf(f,"v %lf %lf %lf\n",vert[i](0),vert[i](1),vert[i](2));
    }
}

void spt::obj_outf(vector<Vector3d> V,vector<Vector3i> F,const char* file_name){
    FILE* f;
    f=fopen(file_name,"w");
    for(int i=0;i<V.size();i++){
        fprintf(f,"v %lf %lf %lf\n",V[i](0),V[i](1),V[i](2));
    }
    for(int i=0;i<F.size();i++){
        fprintf(f,"f %d %d %d\n",F[i](0)+1,F[i](1)+1,F[i](2)+1);
    }
}

// それぞれの柱の高さを設定
void spt::l_height(){
    vector<vector<double>> rays_structure(ohp.size()); // オーバーハング線と交わわるy座標の値
    std::vector<std::vector<Eigen::Vector3d>> rs(ohp.size());
    for(int j=0;j<ohp.size();j++){
        for(int i=0;i<FG.size();i++){
        
            // 内外判定
            Vector3d BA=Vs[FG[i](0)]-Vs[FG[i](1)];
            Vector3d BP=ohp[j]-Vs[FG[i](1)];
            Vector3d BC=Vs[FG[i](2)]-Vs[FG[i](1)];
            Vector3d CP=ohp[j]-Vs[FG[i](2)];
            Vector3d CA=Vs[FG[i](0)]-Vs[FG[i](2)];


            Vector3d cross_a=BP.cross(BA);
            // cout<<"BP*BA={"<<cross_a<<"}"<<endl;
            Vector3d cross_c=BC.cross(BP);
            // cout<<"BC*BP={"<<cross_c<<"}"<<endl;
            Vector3d cross_b=CA.cross(CP);
            // cout<<"CA*CP={"<<cross_b<<"}"<<endl;
            if((cross_a(1)>=0 && cross_b(1)>=0 && cross_c(1)>=0) || (cross_a(1)<=0 && cross_b(1)<=0 && cross_c(1)<=0)){ //内側
            double min=100;
            int ki;
                for(int k=0;k<3;k++){
                    if(VG[FG[i](k)](1)<min){
                        min=VG[FG[i](k)](1);
                        ki=k;
                    }
                }
                rays_structure[j].push_back(min); //oh線と交わった面のy座標
                rs[j].push_back(VG[FG[i](ki)]);
                // cout<<"FG="<<i<<endl;
                // cout<<"ohp="<<ohp[j]<<endl;
                // cout<<"座標:{"<<Vs[FG[i](0)]<<"}\n{"<<Vs[FG[i](1)]<<"}\n{"<<Vs[FG[i](2)]<<"\n"<<endl;
            }
        }
        sort(rays_structure[j].rbegin(),rays_structure[j].rend());
        // for(int i=0;i<rays_structure[j].size();i++){
        //     cout<<"rs="<<rays_structure[j][i]<<endl;
        // }
        height.push_back(rays_structure[j][0]);
    }
    rays_s=rs;
}

// オーバーハング点の判定(入力：グリッドセルの点，座標情報，オーバハング面の構成)
bool spt::point_in_out(Vector3d P,double h){
    for(int i=0;i<oh_fn.size();i++){

        Vector3d BA=Vs[oh_F[i](0)]-Vs[oh_F[i](1)];
        Vector3d BP=P-Vs[oh_F[i](1)];
        Vector3d BC=Vs[oh_F[i](2)]-Vs[oh_F[i](1)];
        Vector3d CP=P-Vs[oh_F[i](2)];
        Vector3d CA=Vs[oh_F[i](0)]-Vs[oh_F[i](2)];


        Vector3d cross_a=BP.cross(BA);
        // cout<<"BP*BA={"<<cross_a<<"}"<<endl;
        Vector3d cross_c=BC.cross(BP);
        // cout<<"BC*BP={"<<cross_c<<"}"<<endl;
        Vector3d cross_b=CA.cross(CP);
        // cout<<"CA*CP={"<<cross_b<<"}"<<endl;

        // cout<<"a="<<cross_a(2)<<" b="<<cross_b(2)<<" c="<<cross_c(2)<<endl;

        if((cross_a(1)>=0 && cross_b(1)>=0 && cross_c(1)>=0) || (cross_a(1)<=0 && cross_b(1)<=0 && cross_c(1)<=0)){ // 外積の向きが揃うとき内側
            // cout<<"内側\n";
            if(VG[oh_F[i](0)](1)-P(1)>h*0.01) return true;
            // 面と地面の距離が閾値以下の場合，オーバーハング点にしない．サポートを立てる必要がないから．
        }

    }
    // cout<<"外側\n";
    return false;
}

// オーバハング面の面番号を設定,元モデルで何番の面なのか
void spt::oh_Fnum(double angle,Vector3d direction){

    for(int i=0;i<FG.size();i++){
        if(overh(angle,direction,VG[FG[i](0)],VG[FG[i](1)],VG[FG[i](2)])){ //オーバハング面だった場合
            oh_fn.push_back(i); //面番号
        }
    }

}

// オーバーハング点の番号を返す(グリッドセルの何番の点か⇒↓)
vector<int> spt::oh_Vnum(vector<Vector3d> gp,double h){
    vector<int> ohv;

    for(int i=0;i<gp.size();i++){
        if(point_in_out(gp[i],h)) ohv.push_back(i); //オーバハング点だった場合
    }
    ohvn=ohv;
    return ohv;
}

//  四隅の点の作成
vector<Vector3d>  spt::f_corners(Vector3d P){
    double x_length=four[1](0)-four[0](0);
    double z_length=four[2](2)-four[0](2);
    double ox=(x_length/N)/2;
    double oz=(z_length/N)/2;

    vector<Vector3d> fc;
    fc.push_back({P(0)-ox,four[0](1),P(2)-oz});
    fc.push_back({P(0)+ox,four[0](1),P(2)-oz});
    fc.push_back({P(0)-ox,four[0](1),P(2)+oz});
    fc.push_back({P(0)+ox,four[0](1),P(2)+oz});

    return fc;
}

// 柱作成(柱の座標，柱の面,底面の頂点座標,柱の高さ,面番号始まりの点)
void spt::pillar(vector<Vector3d> &pv,vector<Vector3i> &pf,vector<Vector3d> fc,double h,int ff){

    for(int i=0;i<4;i++){
        pv.push_back(fc[i]);
    }

    for(int i=0;i<4;i++){
        pv.push_back({fc[i](0),h,fc[i](2)});
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
void spt::L_support(){

    int f_num=0;
    
    for(int i=0;i<ohp.size();i++){
        vector<Vector3d> fc=f_corners(ohp[i]); // 四隅の点を求める

        vector<Vector3d> pv;
        vector<Vector3i> pf;
        // cout<<"height="<<height[i]<<endl;
        pillar(pv,pf,fc,height[i],f_num);
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
    vtom.vmat(VG,FG,v,f);
    cout<<"元モデル変換\n";
    MatrixXd boolsupport_v;
    MatrixXi boolsupport_f;
    igl::MeshBooleanType boolean_type=static_cast<igl::MeshBooleanType>(2);
    cout<<"booleanタイプ選び\n";
    igl::copyleft::cgal::mesh_boolean(spv,spf,v,f,boolean_type,b_v,b_f);
    igl::writeOBJ("boolsupport.obj",b_v,b_f);
}


// int main(int argc,char* argv[])
// {   
//     string fn;
//     FILE *fp;
//     vector<Vector3d> V;
//     vector<Vector3i> F;
//     ro reado;
//     spt sp;

//     if(argc==1){
//         cout<<"入力してください\n";
//         exit(1);
//     }else{
//         fn=argv[1];
//         fp=fopen(argv[1],"r");
//     }
//     reado.readPoint(V,F,fn);
//     sp.VG=V;
//     sp.FG=F;
//     cout<<"読み込み完了\n";
//     // for(int i=0;i<V.size();i++){
//     //     cout<<V[i]<<endl;
//     // }
//     // for(int i=0;i<F.size();i++){
//     //     cout<<F[i]<<endl;
//     // }

//     // x,y平面グリッドセルの作成
//     sp.N=30; //グリッドの分割数
//     sp.aabb(); //最小包含
//     sp.grid_cell();
//     cout<<"グリッドセルの作成完了\n";
//     // obj_out(gc,"grid_cell.obj"); //obj出力

//     double y_max=-100;
//     for(int i=0;i<V.size();i++){ //y座標最大
//         if(y_max<V[i](1)){
//             y_max=V[i](1);
//         }
//     }

//     // 面のオーバハングを調べる
//     cout<<"面のオーバーハングを調べる"<<sp.FG.size()<<endl;
//     double angle=0.5235987755983; //閾値:cura35度(0.6108652381980153)，45度：0.7853981633974483,30度0.5235987755983
//     Vector3d direction(0,-1,0); //造形方向ベクトル
//     vector<int> num; //その面がオーバーハングかどうか(1 or 0)
//     for(int i=0;i<F.size();i++){
//         if(overh(angle,direction,V[F[i](0)],V[F[i](1)],V[F[i](2)])) num.push_back(1);
//         else num.push_back(0);
//     }

//     // オーバハング点の判定
//     sp.oh_Fnum(angle,direction); // オーバーハング面の集合を作る
//     // Vの射影
//     for(int i=0;i<V.size();i++){
//         Vector3d vss={V[i](0),sp.four[0](1),V[i](2)};
//         sp.Vs.push_back(vss);
//     }

//     for(int i=0;i<sp.oh_fn.size();i++){
//         Vector3i ohf={F[sp.oh_fn[i]](0),F[sp.oh_fn[i]](1),F[sp.oh_fn[i]](2)};
//         sp.oh_F.push_back(ohf);
//     }
//     double h=y_max-sp.four[0](1);
//     vector<int> oh_point=sp.oh_Vnum(sp.gc_V,h); //オーバーハング点の点番号
    
//     // デバッグ
//     cout<<"gc="<<sp.gc_V.size()<<" ohp="<<oh_point.size()<<endl;
//     for(int i=0;i<oh_point.size();i++){
//         Vector3d op=sp.gc_V[oh_point[i]];
//         sp.ohp.push_back(op);
//     }
//     cout<<"ohpsize="<<sp.ohp.size()<<endl;

//     sp.obj_out(sp.ohp,"oh_point.obj");

//     vector<Vector3d> sv; //サポート点情報
//     vector<Vector3i> sf; //サポート面情報
//     sp.l_height();
//     cout<<"h="<<h<<endl;
//     sp.L_support(sv,sf);
//     cout<<"sv="<<sv.size()<<" sf="<<sf.size()<<endl;
//     sp.obj_outf(sv,sf,"support_L.obj");

//     // // vtk出力，結果の可視化
//     // FILE* fpv;
//     // fpv=fopen("overhang.vtk","w");
//     // fprintf(fpv,"# vtk DataFile Version 2.0\n");
//     // fprintf(fpv,"Title of my super cool VTK file\nASCII\nDATASET UNSTRUCTURED_GRID\n");
//     // fprintf(fpv,"POINTS %d float\n",(int)V.size());
//     // for(int i=0;i<V.size();i++){
//     //     fprintf(fpv,"%lf %lf %lf\n",V[i](0),V[i](1),V[i](2));
//     // }
//     // fprintf(fpv,"CELLS %d %d\n",(int)F.size(),(int)F.size()*4);
//     // for(int i=0;i<F.size();i++){
//     //     fprintf(fpv,"3 %d %d %d\n",F[i](0),F[i](1),F[i](2));
//     // }
//     // fprintf(fpv,"CELL_TYPES %d\n",(int)F.size());
//     // for(int i=0;i<F.size();i++){
//     //     fprintf(fpv,"5\n");
//     // }
//     // fprintf(fpv,"CELL_DATA %d\n",(int)F.size());
//     // fprintf(fpv,"SCALARS my_face_attr int 1\nLOOKUP_TABLE default\n");
//     // for(int i=0;i<num.size();i++){
//     //     fprintf(fpv,"%d\n",num[i]);
//     // }
//     // cout<<"numsize="<<num.size()<<endl;
// }

// 面のクラス分け
// オーバハングの面だけ残す
// 隣の面がオーバハングか探索する
// https://docs.google.com/spreadsheets/d/1N3K9jsz0g74o5LToadMPtwWfvS1vS7ZjCbrNMR23uYA/edit?usp=sharing
//ビルドする方法：/source$ cmake --build build
//実行する方法：