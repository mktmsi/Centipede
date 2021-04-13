//12-9-1
#include <iostream>
#include <stdio.h>
#include<math.h>
#include <cmath>
#include <stdlib.h>
#include <GLUT/glut.h>
#include <time.h>
#include "Monitor.hpp"
#include "Vector.h"
#include <string.h>
#include <png.h>
#define dt (0.001)	//タイムステップ0.0001
#define MAX_STEP (30000)//最大タイムステップ数
#define PI (3.14159265359)
#define N (30)//質点数

Monitor monitor;

typedef struct{Vector2D r,v;} AA;

AA operator +(AA a,AA b){AA c;	c.r=a.r+b.r;	c.v=a.v+b.v;	return c;}
AA operator -(AA a,AA b){AA c;	c.r=a.r-b.r;	c.v=a.v-b.v ;	return c;}
AA ax(double a,AA b){AA c;	c.r=b.r*a;		c.v=b.v*a;	return c;}
AA wx(double a,AA b){AA c;	c.r=b.r/a;		c.v=b.v/a;	return c;}

#define aa(i,j) aa[N*(j)+(i)]
#define r(i,j) r[N*(j)+(i)]
#define v(i,j) v[N*(j)+(i)]

AA aa[6*N];

//質点の位置・相互関係
Vector2D l[N];//質点間の方向ベクトル
Vector2D et[N];//質点間の単位方向ベクトル
Vector2D en[N];//質点間の単位法線ベクトル
Vector2D et_hat[N];//体軸方向の単位ベクトル
Vector2D en_hat[N];//体軸方向に対して垂直な単位方向ベクトル
double absl[N];//質点間の距離
double ldot[N];//質点間距離の変化量
double lbar=2.0;
double l0=2.0;
double phi[N],phidot[N],phibar[N];//i番目の質点の角度

//質点にかかる力・トルク
Vector2D Fbody[N];//バネ・ダンパ由来の力ベクトル
Vector2D Ftorque[N];//トルク由来の力ベクトル
Vector2D Ffriction[N];//摩擦抵抗の力ベクトル
double fs[N],fd[N];
double tau[N],tau_pas[N],tau_act[N];

//機械定数
double k=1000.0,kpas=5.1,kact=20.0;//バネ定数kpas=15.1,kact=50.0
double c=100.5,cpas=1.11;//ダンパ定数cpas=1.1
double mt=0.001,mn=1.0;//摩擦定数 mn=0.001
double m=5.0;//質量5.0

//先頭関節のトルクa*sin(wt)
double a=1.0;//a=2.0
double w=20.0;//100.0;

double t;
double L=100.0;
int winid;
int ts;//現在のタイムステップ
int initial=1;
FILE *fp;	//fpというファイル用変数を定義


void capture(int *);

const int save_flag=0;    //0:画像保存しない，1:画像保存する

void init(){
  int i,j;
  double theta;
  char filename[256];
  t,ts=0;
  //等方性摩擦
  //mt=mn;
  sprintf(filename, "data.txt");    //数値データ保存ファイルのファイル名設定
  fp=fopen(filename,"w");        //数値データ保存ファイルを開く

  printf("initial place was initialized\n");
  for(i=0;i<N;i++){
    aa(i,0).r.set_x(-l0*i);
    aa(i,0).r.set_y(0.0);//r(i,0).set_x();
    //r(i,0).set_x(-l0*i);
    //r(i,0).set_y(0.0);//r(i,0).set_x();
  }
}


void func(AA *p, AA *p_out){
  int i;
  Vector2D temp;
  double kaku;

  //各質点にかかるバネ・ダンパ由来の力ベクトルを更新
  for(i=1;i<N;i++){
    fs[i] = k*pow(absl[i]-lbar,3);//spring
    fd[i] = c*ldot[i];//damper
  }
  for(i=0;i<N;i++){
    if(i==0){
      Fbody[i] = et[i+1]*(-(fs[i+1]+fd[i+1]));
    }else if(i==N-1){
      Fbody[i] = et[i]*(fs[i]+fd[i]);
    }else{
      Fbody[i] = et[i]*(fs[i]+fd[i])+et[i+1]*(-(fs[i+1]+fd[i+1]));
    }
  }

  //質点の角度を更新
  for(i=1;i<N;i++){
    if(i==N-1){
      phi[i] = atan2(p[i-1].r.get_y()-p[i].r.get_y(),p[i-1].r.get_x()-p[i].r.get_x());
      phidot[i] = (p[i-1].v-p[i].v)*en[i]/absl[i];
    }else{
      phi[i] = atan2(et[i+1]^et[i],et[i+1]*et[i]);
      phidot[i] = ((p[i-1].v-p[i].v)*en[i])/absl[i] + ((p[i+1].v-p[i].v)*en[i+1])/absl[i+1];
    }
  }

  //目標角度・質点にかかるトルクを更新
  for(i=1;i<N-1;i++){
    if(i==1){
      tau[i] = a*sin(w*t);
      tau[i] += -kpas*phi[i] - cpas*phidot[i];
    }else{
      phibar[i] = phi[i-1];
      tau_pas[i] = -kpas*phi[i] - cpas*phidot[i];
      tau_act[i] = -kact*(phi[i] - phibar[i]);
      tau[i] = tau_pas[i] + tau_act[i];
    }
  }

  //質点にかかるトルク由来の力ベクトルを更新
  for(i=0;i<N;i++){
    if(i==0){
      Ftorque[0] = en[1]*(tau[1]/absl[1]);
    }else if(i==1){
      Ftorque[1] = en[1]*(-tau[1]/absl[1]) + en[2]*((tau[2]-tau[1])/absl[2]);
    }else if(i==N-2){
      Ftorque[N-2] = en[N-2]*((tau[N-3]-tau[N-2])/absl[N-2]) + en[N-1]*(-tau[N-2]/absl[N-1]);
    }else if(i==N-1){
      Ftorque[N-1] = en[N-1]*(tau[N-2]/absl[N-1]);
    }else{
      Ftorque[i] = en[i]*((tau[i-1]-tau[i])/absl[i]) + en[i+1]*((tau[i+1]-tau[i])/absl[i+1]);
    }
  }

  //摩擦力を更新
  for(i=0;i<N;i++){
    Ffriction[i] = et_hat[i]*((-mt)*(p[i].v*et_hat[i])) + en_hat[i]*((-mn)*(p[i].v*en_hat[i]));
  }

  //運動方程式
  for(i=0;i<N;i++){
    p_out[i].v=(Fbody[i] + Ftorque[i] + Ffriction[i])*(dt/m);
    p_out[i].r = p[i].v*dt;
  }

}

void UpdateUnitVectors(AA *p){
  int i;
  Vector2D temp;

  //質点の位置関係に関する変数を更新
  for(i=1;i<N;i++){
    l[i] = p[i-1].r - p[i].r;//質点間の相対距離
    absl[i] = l[i].get_abs();
    et[i] = l[i] / absl[i];
    en[i].set_x(et[i].get_y()*sin(-PI/2));//回転行列
    en[i].set_y(et[i].get_x()*sin(PI/2));//回転行列
    ldot[i] = (p[i-1].v-p[i].v) * et[i];//体軸方向の相対速度性分
  }

  //体軸方向の単位ベクトルを更新
  for(i=1;i<N-1;i++){
    temp = et[i] + et[i+1];
    et_hat[i] = temp/temp.get_abs();
    en_hat[i].set_x(et_hat[i].get_y()*sin(-PI/2));//回転行列
    en_hat[i].set_y(et_hat[i].get_x()*sin(PI/2));//回転行列
  }
  //i=0のときの体軸方向の単位ベクトルを更新
  et_hat[0].set_x(et_hat[1].get_x()*cos(phi[1])+et_hat[1].get_y()*sin(-phi[1]));//回転行列
  et_hat[0].set_y(et_hat[1].get_x()*sin(phi[1])+et_hat[1].get_y()*cos(phi[1]));//回転行列
  en_hat[0].set_x(et_hat[0].get_y()*sin(-PI/2));//回転行列
  en_hat[0].set_y(et_hat[0].get_x()*sin(PI/2));//回転行列

  //i=N-1のときの体軸方向の単位ベクトルを更新
  et_hat[N-1].set_x(et_hat[N-2].get_x()*cos(-phi[N-2])+et_hat[N-2].get_y()*sin(phi[N-2]));//回転行列
  et_hat[N-1].set_y(et_hat[N-2].get_x()*sin(-phi[N-2])+et_hat[N-2].get_y()*cos(-phi[N-2]));//回転行列
  en_hat[N-1].set_x(et_hat[N-1].get_y()*sin(-PI/2));//回転行列
  en_hat[N-1].set_y(et_hat[N-1].get_x()*sin(PI/2));//回転行列
}

void runge(){
  int i,j,q,n;
  //i番目の質点の位置に関する微分方程式を解く
  for (q=0; q<100; q++) {
    UpdateUnitVectors(&aa(0,0));
    func(&aa(0,0),&aa(0,1));    //aa(0,1)=k1
    for (j=0; j<N; j++) {
      aa(j,5)=aa(j,0)+wx(2.0,aa(j,1));
    }
    UpdateUnitVectors(&aa(0,5));
    func(&aa(0,5),&aa(0,2));    //aa(0,2)=k2
    for (j=0; j<N; j++) {
      aa(j,5)=aa(j,0)+wx(2.0,aa(j,2));
    }
    UpdateUnitVectors(&aa(0,5));
    func(&aa(0,5),&aa(0,3));    //aa(0,3)=k3
    for (j=0; j<N; j++) {
      aa(j,5)=aa(j,0)+aa(j,3);
    }
    UpdateUnitVectors(&aa(0,5));
    func(&aa(0,5),&aa(0,4));    //aa(0,4)=k4
    for (j=0; j<N; j++) {
      aa(j,0)=aa(j,0)+wx(6.0,aa(j,1)+ax(2.0,aa(j,2))+ax(2.0,aa(j,3))+aa(j,4)); //x(t+dt)=x(t)+k
    }
  }

  ts++;
  t=ts*dt;
  //結果の表示
  if(ts%100==0){
    printf("t=%f  \n",(double)t);
    if(ts>= MAX_STEP){ //MAX_STEPに達したらプログラムを終了
      fclose(fp);
      exit(0);
    }
  }
}



void keyboard(unsigned char key, int x , int y){
  double tmpcin;
  switch(key)
  {
    case 'Q'   :                                break;
    case 'q'   :                                break;
    case '\033':  /* '\033' = ESC */ exit(0);   break;
    case 'h'   : monitor.SetCenter( 1.0 , 0 );  break;
    case 'l'   : monitor.SetCenter( -1.0 , 0 ); break;
    case 'j'   : monitor.SetCenter( 0 , 1.0 );  break;
    case 'k'   : monitor.SetCenter( 0 , -1.0 ); break;
    case 'c'   :                                break;
    case 'g'   :                                break;
    case 'z'   : monitor.SetZoom( 1.1/ 1.0 );   break;
    case 'x'   : monitor.SetZoom( 1.0/ 1.1 );  break;
  }

}

void idle(void)
{
  glutSetWindow(winid);
  //glutKeyboardFunc(keyboard);
  glutPostRedisplay();
}

void display(void)
{
  glClear ( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
  int i,j;
  double o=0.0,n;
  //ディスプレイ表示用文字列
  //char str[256],str1[256],str2[256],str3[256],str4[256],str5[256];
  char str[256];
  glClear(GL_COLOR_BUFFER_BIT);
  glLoadIdentity();
  //glutKeyboardFunc(keyboard);

  runge();        //runge-kuttaを回す

  /*ィンドウサイズに表示*/
  monitor.SetCenter(aa(j,0).r.get_x()/L*2-1.5  , 0 );
  //  o=n;
  for(j=0;j<N;j++){
    if(tau[j]>0){//トルクが正であればその質点は青色
      monitor.SetAllColor(0.0,0.0,1.0);//青
    }else{
      monitor.SetAllColor(1.0,0.0,.00);//赤
    }
    monitor.DrawCircle(aa(j,0).r.get_x()/L*2-1.5,aa(j,0).r.get_y()/L*2,0.02);
  }

  monitor.SetAllColor(0.0,0.0,0.0);
  for(j=0;j<N;j++){
    if(j!=N-1){
      monitor.DrawLine(aa(j,0).r.get_x()/L*2-1.5,aa(j,0).r.get_y()/L*2,aa(j+1,0).r.get_x()/L*2-1.5,aa(j+1,0).r.get_y()/L*2,0.005);

    }
  }
  for(j=-N/3;j<300;j++){
    monitor.DrawLine(0.3*j,1,0.3*(j),-2,0.005);
  }
  monitor.SetAllColor(0.0,0.0,0.0);

  sprintf(str,"Isotropic friction");
  monitor.String(0.1,0.6,str);
  glFlush();
  glutSwapBuffers();

  //1000ステップごとに画像を保存
  if(!(ts%100)) {
    capture(&ts);
  }
}




void mouse(int button, int state, int x, int y){
  switch (button)
  {
    case GLUT_LEFT_BUTTON:
    if (state == GLUT_DOWN) { glutIdleFunc(0);     std::cout << "left: on" << std::endl;  }
    else                    { glutIdleFunc(idle);  std::cout << "left: off" << std::endl; }
    break;

    case GLUT_MIDDLE_BUTTON:
    if (state == GLUT_DOWN) { glutIdleFunc(0);    std::cout << "middle: on" << std::endl;  }
    else                    { glutIdleFunc(idle); std::cout << "middle: off" << std::endl; }
    break;

    case GLUT_RIGHT_BUTTON:
    if (state == GLUT_DOWN) { glutIdleFunc(0);    std::cout << "right: on" << std::endl;  }
    else                    { glutIdleFunc(idle); std::cout << "right: off" << std::endl; }
    break;
  }
}

void resize( int w , int h )
{
  glViewport(0, 0, w, h);

  //glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  glOrtho( -w / 0.0 , w / 20.0 , -h / 20.0 , h / 20.0 , -3.0 , 3.0 );
  //gluPerspective( 30.0, (double)w / (double)h, 1.0 , 100.0 );
  //gluLookAt( 0.0 , 0.0 , 3.8 , 0.0 , 0.0 , 0.0 , 0.0 , 1.0 , 0.0 );
  //monitor.SetWindowSize( w , h );
  //glMatrixMode(GL_MODELVIEW);
}



void OpenGL_init(int *argcp , char **argv)
{
  init();        //初期条件を設定

  glutInit(argcp, argv);

  glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE);

  glutInitWindowSize(monitor.GetWindowSize(Monitor::X),monitor.GetWindowSize(Monitor::Y));
  glutInitWindowPosition( 10 , 100 );//(10,100)
  winid = glutCreateWindow("simulation");
  glutDisplayFunc(display);
  glutReshapeFunc(resize);
  glutKeyboardFunc(keyboard);
  glutMouseFunc(mouse);
  glClearColor(1.0, 1.0, 1.0, 1.0);
}

void monitor_init()
{
  monitor.SetWindowSize( 800 , 600 );
  //   monitor.SetMode( 0 );
  monitor.SetMovieMode( 1);
  monitor.SetMovieName( "./MovieDir/temp_" );
  //   monitor.SetGridMode( 0 );
  //   monitor.SetGridWidth( 2.0 );

}

void capture(int *pts)
{
  char filepath[100]; //= "./MovieDir/output.png";
  sprintf(filepath, "./MovieDir/%d.png", *pts);
  png_bytep raw1D;
  png_bytepp raw2D;
  int i;
  int width = glutGet(GLUT_WINDOW_WIDTH);
  int height = glutGet(GLUT_WINDOW_HEIGHT);

  // 構造体確保
  FILE *fp = fopen(filepath, "wb");
  png_structp pp = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
  png_infop ip = png_create_info_struct(pp);
  // 書き込み準備
  png_init_io(pp, fp);
  png_set_IHDR(pp, ip, width, height,
               8,                   // 8bit以外にするなら変える
               PNG_COLOR_TYPE_RGBA, // RGBA以外にするなら変える
               PNG_INTERLACE_NONE, PNG_COMPRESSION_TYPE_DEFAULT, PNG_FILTER_TYPE_DEFAULT);
  // ピクセル領域確保
  raw1D = (png_bytep)malloc(height * png_get_rowbytes(pp, ip));
  raw2D = (png_bytepp)malloc(height * sizeof(png_bytep));
  for (i = 0; i < height; i++)
    raw2D[i] = &raw1D[i * png_get_rowbytes(pp, ip)];
  // 画像のキャプチャ
  glReadBuffer(GL_FRONT);
  glPixelStorei(GL_UNPACK_ALIGNMENT, 1); // 初期値は4
  glReadPixels(0, 0, width, height,
               GL_RGBA,          // RGBA以外にするなら変える
               GL_UNSIGNED_BYTE, // 8bit以外にするなら変える
               (void *)raw1D);
  // 上下反転
  for (i = 0; i < height / 2; i++)
  {
    png_bytep swp = raw2D[i];
    raw2D[i] = raw2D[height - i - 1];
    raw2D[height - i - 1] = swp;
  }
  // 書き込み
  png_write_info(pp, ip);
  png_write_image(pp, raw2D);
  png_write_end(pp, ip);
  // 開放
  png_destroy_write_struct(&pp, &ip);
  fclose(fp);
  free(raw1D);
  free(raw2D);

  //printf("write out screen capture to '%s'\n", filepath);
}


int main(int argc, char *argv[])
{
  int i, j, i_dim;

  monitor_init();
  std::cout << "monitor init OK" << std::endl;

  OpenGL_init(&argc, argv);

  std::cout << "OpenGL init OK" << std::endl;

  // glutKeyboardFunc(keyboard);
  glutIdleFunc(idle);
  glutMainLoop();//無限ループ

  return 0;
}
