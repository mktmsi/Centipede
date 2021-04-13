//12-9-1
#include <iostream>
#include <stdio.h>
#include <math.h>
#include <cmath>
#include <stdlib.h>
#include <GLUT/glut.h>
#include <time.h>
#include "Monitor.hpp"
#include "Vector.h"
#include <string.h>
#include <png.h>
#define dt (0.0001)      //タイムステップ0.0001
#define MAX_STEP (30000) //最大タイムステップ数
#define PI (3.14159265359)
#define N (10) //質点数
#define trunk(i, j) trunk[N * (j) + (i)]
#define leg(i, j) leg[N * (j) + (i)]
#define phi_osci(i, j) phi_osci[N * (j) + (i)]
//#define r(i, j) r[N * (j) + (i)]
//#define v(i, j) v[N * (j) + (i)]

Monitor monitor;

typedef struct
{
  Vector2D r, v;
} AA;
AA operator+(AA a, AA b)
{
  AA c;
  c.r = a.r + b.r;
  c.v = a.v + b.v;
  return c;
}
AA operator-(AA a, AA b)
{
  AA c;
  c.r = a.r - b.r;
  c.v = a.v - b.v;
  return c;
}
AA ax(double a, AA b)
{
  AA c;
  c.r = b.r * a;
  c.v = b.v * a;
  return c;
}
AA wx(double a, AA b)
{
  AA c;
  c.r = b.r / a;
  c.v = b.v / a;
  return c;
}

//AA trunk[6*N];

/*********************質点********************/
class Mass
{
private:
public:
  double m, battery;
  //Vector2D F;
  AA aa;
  //int state_num;
};

Mass trunk[6 * (N + 1)];
Mass leg[6 * (N + 1)];

//質点の位置・相互関係
Vector2D vecttor_tt[N], vecttor_tl[N];     //質点間の方向ベクトル
Vector2D et_tt[N], et_tl[N];               //質点間の単位方向ベクトル
Vector2D en_tt[N], en_tl[N];               //質点間の単位法線ベクトル
Vector2D et_tt_hat[N], et_tl_hat[N];       //体軸方向の単位ベクトル
Vector2D en_tt_hat[N], en_tl_hat[N];       //体軸方向に対して垂直な単位方向ベクトル
double abs_tt[N], abs_tl[N];               //質点間の距離
double l_ttDot[N], l_tlDot[N];             //質点間距離の変化量
double l_tbar = 10.0, l_lbar = 10.0;       //体節間の自然長
double l_t0 = l_tbar, l_l0 = l_lbar;       //体節間の初期長，脚の書記長
double phi_t[N], phi_tDot[N], phi_tBar[N]; //i番目の質点の角度
double phi_l[N], phi_lDot[N], phi_lBar[N]; //i番目の質点の角度

//オシレックス
double phi_osci[N * 6], lbar_osci[N], l0_osci = l_l0;
double Bstan = 2.0, Bswin = 008.0, alpha = 1.0;
const double omega = 2.0, sigma = 0.15 * 3.0;

//質点にかかる力・トルク
Vector2D F_tBody[N], F_lBody[N];         //バネ・ダンパ由来の力ベクトル
Vector2D F_tTorque[N], F_lTorque[N];     //トルク由来の力ベクトル
Vector2D F_tFriction[N], F_lFriction[N]; //摩擦抵抗の力ベクトル
Vector2D F_lGround[N], F_tGround[N];
double ft_s[N], ft_d[N], fl_s[N], fl_d[N];
double tau_t[N], tau_tPas[N], tau_tAct[N], tau_l[N], tau_lPas[N], tau_lAct[N];

//機械定数
double kt = 10000.0, ct = 100.5, kl = 1000.0, cl = 100000.5; //直動バネ・ダンパ                      //
double kt_pas = 10000.1, kt_act = 100000.0;                  //まきダンパ　バネ定数kpas=15.1,kt_act=50.0
double kl_pas = 1000.1, kl_act = 10000.0;                    //double kl_pas = 1000.1, kl_act = 2000000.0;
double ct_pas = 10000.11;                                    //ダンパ定数cpas=1.1
double cl_pas = 10000.11;                                    //ダンパ定数cpas=1.1
//double mt = 0.001, mn = 1.0;      //摩擦定数 mn=0.001
double m = 10.0; //質量5.0
const double c_gr = 100.0, k_gr = 100000.0;
const double mu = 0.5 * 1.0;
const double c_tanh = 10.0;

double t;
double L = 100.0;
int winid;
int ts; //現在のタイムステップ
int initial = 1;
FILE *fp; //fpというファイル用変数を定義

void capture(int *);

const int save_flag = 0; //0:画像保存しない，1:画像保存する

void init()
{
  int i, j;
  double theta;
  char filename[256];
  t, ts = 0;
  //等方性摩擦
  //mt=mn;
  sprintf(filename, "data.txt"); //数値データ保存ファイルのファイル名設定
  fp = fopen(filename, "w");     //数値データ保存ファイルを開く

  printf("initial place was initialized\n");
  for (i = 0; i < N; i++)
  {
    leg(i, 0).aa.r.set_x(-l_t0 * i + 0.0);
    leg(i, 0).aa.r.set_y(0.0); //r(i,0).set_x();
    trunk(i, 0).aa.r.set_x(-l_t0 * i);
    trunk(i, 0).aa.r.set_y(leg(i, 0).aa.r.get_y() + l_l0); //r(i,0).set_x();
    //trunk(i, 0).m=m;
    //leg(i, 0).m=m;
    trunk(i, 0).aa.v.set_vec(0.0, 0.0);
    leg(i, 0).aa.v.set_vec(0.0, 0.0); //r(i,0).set_x();
    phi_osci(i, 0) = PI / 2.0;        //r(i,0).set_x();
    lbar_osci[i] = l0_osci - Bswin * max(0.0, sin(phi_osci(i, 0))) - Bstan * min(0.0, sin(phi_osci(i, 0)));
    phi_lBar[i] = -alpha * cos(phi_osci(i, 0)); //phi_l[i - 1];
  }
  /*trunk(0, 0).m=0.0;
  trunk(N-1, 0).m=0.0;
  leg(0, 0).m=0.0;
  leg(N-1, 0).m=0.0;
  */

  printf("x0_x=%lf,v0_x=%lf\n", trunk[0].aa.r.get_x(), trunk[0].aa.v.get_x());
}

void func(Mass *pt, Mass *pt_out, Mass *pl, Mass *pl_out, double *pPhi, double *pPhi_out)
{
  int i;
  Vector2D temp_t, temp_l;
  Vector2D g;
  g.set_vec(0.0, -9.8);
  double kaku;

  //足の長さを更新，目標角度を更新
  for (i = 1; i < N; i++)
  {
    lbar_osci[i] = l0_osci - Bswin * max(0.0, sin(pPhi[i])) - Bstan * min(0.0, sin(pPhi[i]));
    phi_lBar[i] = -alpha * cos(pPhi[i]); //phi_l[i - 1];
  }

  //質点の位置関係に関する変数を更新
  for (i = 1; i < N; i++)
  {
    //体節間の相対関係
    vecttor_tt[i] = pt[i - 1].aa.r - pt[i].aa.r; //質点間の相対距離
    abs_tt[i] = vecttor_tt[i].get_abs();
    et_tt[i] = vecttor_tt[i] / abs_tt[i];
    en_tt[i].set_x(et_tt[i].get_y() * sin(-PI / 2));       //回転行列
    en_tt[i].set_y(et_tt[i].get_x() * sin(PI / 2));        //回転行列
    l_ttDot[i] = (pt[i - 1].aa.v - pt[i].aa.v) * et_tt[i]; //体軸方向の相対速度性分
    //体節-脚間の相対関係

    vecttor_tl[i] = pl[i].aa.r - pt[i].aa.r; //質点間の相対距離
    abs_tl[i] = vecttor_tl[i].get_abs();
    et_tl[i] = vecttor_tl[i] / abs_tl[i];
    en_tl[i].set_x(et_tl[i].get_y() * sin(-PI / 2));   //回転行列
    en_tl[i].set_y(et_tl[i].get_x() * sin(PI / 2));    //回転行列
    l_tlDot[i] = (pl[i].aa.v - pt[i].aa.v) * et_tl[i]; //体軸方向の相対速度性分
  }

  //各質点にかかるバネ・ダンパ由来の力ベクトルを更新
  for (i = 0; i < N; i++)
  {
    if (i == 0) //あってる
    {
      ft_s[i] = 0.0;
      ft_d[i] = 0.0;
    }
    else
    {
      ft_s[i] = kt * pow(abs_tt[i] - l_tbar, 3); //spring
      ft_d[i] = ct * l_ttDot[i];                 //damper
    }

    if (i == 0 || i == N - 1)
    {
      fl_s[i] = 0.0;
      fl_d[i] = 0.0;
    }
    else
    {
      fl_s[i] = kl * pow(abs_tl[i] - lbar_osci[i], 3); //spring
      fl_d[i] = cl * l_tlDot[i];                       //damper
    }
  }

  for (i = 0; i < N; i++)
  {
    if (i == 0)
    {
      F_tBody[i] = et_tt[i + 1] * (-(ft_s[i + 1] + ft_d[i + 1]));
      F_lBody[i].set_vec(0.0, 0.0);
      //F_lBody[i] = et_tl[i+1]*(-(fl_s[i+1]+fl_d[i+1]));
    }
    else if (i == N - 1)
    {
      F_tBody[i] = et_tt[i] * (ft_s[i] + ft_d[i]);
      F_lBody[i].set_vec(0.0, 0.0);
      //F_lBody[i] = et_tl[i]*(fl_s[i]+fl_d[i]);
    }
    else
    {
      F_tBody[i] = et_tt[i] * (ft_s[i] + ft_d[i]) + et_tt[i + 1] * (-(ft_s[i + 1] + ft_d[i + 1])) + et_tl[i] * (fl_s[i] + fl_d[i]);
      F_lBody[i] = (et_tl[i]) * (-(fl_s[i] + fl_d[i]));
    }
  }

  //質点の角度を更新
  for (i = 1; i < N; i++)
  {
    //体節
    if (i == N - 1)
    {
      phi_t[i] = atan2(pt[i - 1].aa.r.get_y() - pt[i].aa.r.get_y(), pt[i - 1].aa.r.get_x() - pt[i].aa.r.get_x());
      phi_tDot[i] = (pt[i - 1].aa.v - pt[i].aa.v) * en_tt[i] / abs_tt[i];
    }
    else
    {
      phi_t[i] = atan2(et_tt[i + 1] ^ et_tt[i], et_tt[i + 1] * et_tt[i]);
      phi_tDot[i] = ((pt[i - 1].aa.v - pt[i].aa.v) * en_tt[i]) / abs_tt[i] + ((pt[i + 1].aa.v - pt[i].aa.v) * en_tt[i + 1]) / abs_tt[i + 1];
    }
    /*if (phi_t[i] > 2.0 * PI)
    {
      phi_t[i] -= 2.0 * PI;
    }
    else if (phi_t[i] < 0.0 )
    {
      phi_t[i] += 2.0 * PI;
    }*/
    //脚
    phi_l[i] = atan2(et_tt[i + 1] ^ et_tl[i], et_tt[i + 1] * et_tl[i]) - PI / 2.0 - PI; //ホントはpi/2だけ引けばいいけど，プログラムの仕様のせいか-PIもしないといけない
    if (phi_l[i] > 2.0 * PI)
    {
      phi_l[i] -= 2.0 * PI;
    }
    if (phi_l[i] < 0.0)
    {
      phi_l[i] += 2.0 * PI;
    }
    phi_lDot[i] = ((pl[i].aa.v - pt[i].aa.v) * en_tl[i]) / abs_tl[i] + ((pt[i + 1].aa.v - pt[i].aa.v) * en_tt[i + 1]) / abs_tt[i + 1];
  }

  //目標角度・質点にかかるトルクを更新
  for (i = 1; i < N - 1; i++)
  {
    //体節
    phi_tBar[i] = 0.0;
    tau_tPas[i] = -kt_pas * phi_t[i] - ct_pas * phi_tDot[i];
    tau_tAct[i] = -kt_act * (phi_t[i] - phi_tBar[i]);
    tau_t[i] = tau_tPas[i] + tau_tAct[i];

    //脚
    tau_lPas[i] = -kl_pas * phi_l[i] - cl_pas * phi_lDot[i];
    tau_lAct[i] = -kl_act * (phi_l[i] - phi_lBar[i]);
    tau_l[i] = tau_lPas[i] + tau_lAct[i];
  }

  //質点にかかるトルク由来の力ベクトルを更新
  for (i = 0; i < N; i++)
  {

    //体節
    if (i == 0)
    {
      F_tTorque[0] = en_tt[1] * (tau_t[1] / abs_tt[1]);
    }
    else if (i == 1)
    {
      F_tTorque[1] = en_tt[1] * (-tau_t[1] / abs_tt[1]) + en_tt[2] * ((tau_t[2] - tau_t[1]) / abs_tt[2]);
    }
    else if (i == N - 2)
    {
      F_tTorque[N - 2] = en_tt[N - 2] * ((tau_t[N - 3] - tau_t[N - 2]) / abs_tt[N - 2]) + en_tt[N - 1] * (-tau_t[N - 2] / abs_tt[N - 1]);
    }
    else if (i == N - 1)
    {
      F_tTorque[N - 1] = en_tt[N - 1] * (tau_t[N - 2] / abs_tt[N - 1]);
    }
    else
    {
      F_tTorque[i] = en_tt[i] * ((tau_t[i - 1] - tau_t[i]) / abs_tt[i]) + en_tt[i + 1] * ((tau_t[i + 1] - tau_t[i]) / abs_tt[i + 1]);
    }
    //脚
    if (i == 0 || i == N - 1)
    {
      F_lTorque[i].set_vec(0.0, 0.0);
    }
    else
    {
      F_lTorque[i] = (en_tl[i] * tau_l[i]) / abs_tl[i];
    }
    F_tTorque[i] -= F_lTorque[i];
  }

  double xabs = 0.0, tanX = 0.0;
  //床からの粘弾性反力と摩擦力
  for (i = 0; i < N - 1; i++)
  {
    if (trunk[i].aa.r.get_y() < 0.0)
    {
      F_tGround[i].set_x(0.0);
      F_tGround[i].set_y(max(-k_gr * trunk[i].aa.r.get_y() - c_gr * trunk[i].aa.v.get_y(), 0.0));
    }
    else
    {
      F_tGround[i].set_vec(0.0, 0.0);
    }
    if (leg[i].aa.r.get_y() < 0.0)
    {
      F_lGround[i].set_x(0.0);
      F_lGround[i].set_y(max(-k_gr * leg[i].aa.r.get_y() - c_gr * leg[i].aa.v.get_y(), 0.0));
    }
    else
    {
      F_lGround[i].set_vec(0.0, 0.0);
    }

    xabs = fabs(F_lGround[i].get_y());
    tanX = tanh(c_tanh * leg[i].aa.v.get_x());
    F_lFriction[i].set_vec(-mu * xabs * tanX, 0.0);
    xabs = fabs(F_tGround[i].get_y());
    tanX = tanh(c_tanh * trunk[i].aa.v.get_x());
    F_tFriction[i].set_vec(-mu * xabs * tanX, 0.0);
  }

  F_lGround[0].set_vec(0.0, 0.0);
  F_lFriction[0].set_vec(0.0, 0.0);
  F_lGround[N - 1].set_vec(0.0, 0.0);
  F_lFriction[N - 1].set_vec(0.0, 0.0);

  /*位相時間発展式*/
  for (i = 1; i < N - 1; i++)
  {
    //PHI_OUT[i]=dt*(omega-sigma_v*F_grf_v[4].y*cos(PHI[i])+sigma_spine*max(0,Ts[2]+Td[2])*cos(PHI[i]));
    pPhi_out[i] = dt * (omega - sigma * tanh(F_lGround[i].get_y()) * cos(pPhi[i]));
  }
  //運動方程式
  for (i = 0; i < N; i++)
  {
    pt_out[i].aa.v = (F_tBody[i] + F_tTorque[i] + F_tFriction[i] + F_tGround[i] + g * m) * (dt / m);
    pt_out[i].aa.r = pt[i].aa.v * dt;
    pl_out[i].aa.v = (F_lBody[i] + F_lTorque[i] + F_lFriction[i] + F_lGround[i] + g * m) * (dt / m);
    pl_out[i].aa.r = pl[i].aa.v * dt;
  }
}

void runge()
{
  int i, j, q, n;
  //i番目の質点の位置に関する微分方程式を解く
  for (q = 0; q < 100; q++)
  {
    //printf("runge  x0_x=%lf,v0_x=%lf\n",trunk(0,0).aa.r.get_x(),trunk(0,0).aa.v.get_x());
    //UpdateUnitVectors(&trunk(0,0).aa);
    //printf("runge2  x0_x=%lf,v0_x=%lf\n",trunk(0,0).aa.r.get_x(),trunk(0,0).aa.v.get_x());
    func(&trunk(0, 0), &trunk(0, 1), &leg(0, 0), &leg(0, 1), &phi_osci(0, 0), &phi_osci(0, 1)); //trunk(0,1)=k1
    for (j = 0; j < N; j++)
    {
      trunk(j, 5).aa = trunk(j, 0).aa + wx(2.0, trunk(j, 1).aa);
      leg(j, 5).aa = leg(j, 0).aa + wx(2.0, leg(j, 1).aa);
      phi_osci(j, 5) = phi_osci(j, 0) + phi_osci(j, 1) / 2.0;
    }
    //UpdateUnitVectors(&trunk(0,5).aa);
    func(&trunk(0, 5), &trunk(0, 2), &leg(0, 5), &leg(0, 2), &phi_osci(0, 5), &phi_osci(0, 2)); //trunk(0,2)=k2
    for (j = 0; j < N; j++)
    {
      trunk(j, 5).aa = trunk(j, 0).aa + wx(2.0, trunk(j, 2).aa);
      leg(j, 5).aa = leg(j, 0).aa + wx(2.0, leg(j, 2).aa);
      phi_osci(j, 5) = phi_osci(j, 0) + phi_osci(j, 2) / 2.0;
    }
    //UpdateUnitVectors(&trunk(0,5).aa);
    func(&trunk(0, 5), &trunk(0, 3), &leg(0, 5), &leg(0, 3), &phi_osci(0, 5), &phi_osci(0, 3)); //trunk(0,3)=k3
    for (j = 0; j < N; j++)
    {
      trunk(j, 5).aa = trunk(j, 0).aa + trunk(j, 3).aa;
      leg(j, 5).aa = leg(j, 0).aa + leg(j, 3).aa;
      phi_osci(j, 5) = phi_osci(j, 0) + phi_osci(j, 3);
    }
    //UpdateUnitVectors(&trunk(0,5).aa);
    func(&trunk(0, 5), &trunk(0, 4), &leg(0, 5), &leg(0, 4), &phi_osci(0, 5), &phi_osci(0, 4)); //trunk(0,4)=k4
    for (j = 0; j < N; j++)
    {
      trunk(j, 0).aa = trunk(j, 0).aa + wx(6.0, trunk(j, 1).aa + ax(2.0, trunk(j, 2).aa) + ax(2.0, trunk(j, 3).aa) + trunk(j, 4).aa); //x(t+dt)=x(t)+k
      leg(j, 0).aa = leg(j, 0).aa + wx(6.0, leg(j, 1).aa + ax(2.0, leg(j, 2).aa) + ax(2.0, leg(j, 3).aa) + leg(j, 4).aa);             //x(t+dt)=x(t)+k
      phi_osci(j, 0) = phi_osci(j, 0) + (phi_osci(j, 1) + 2.0 * phi_osci(j, 2) + 2.0 * phi_osci(j, 3) + phi_osci(j, 4)) / 6.0;
      if (phi_osci(j, 0) > 2.0 * PI)
      {
        phi_osci(j, 0) -= 2.0 * PI;
      }
      if (phi_osci(j, 0) < 0.0)
      {
        phi_osci(j, 0) += 2.0 * PI;
      }
    }
  }
  /*printf("ts=%d\n",ts);
  printf("runge  x0_lx=%lf,v0_lx=%lf  y0_lx=%lf,v0_ly=%lf\n", leg(1, 0).aa.r.get_x(), leg(1, 0).aa.v.get_x(), leg(1, 0).aa.r.get_y(), leg(1, 0).aa.v.get_y());
  printf("runge  phi_t1=%lf,  phi_l1=%lf\n", phi_t[1],phi_l[1]);
  printf("runge  Torque0_tx=%lf,  Torque0_ty=%lf\n", F_tTorque[0].get_x(),F_tTorque[0].get_y());*/
  ts++;
  t = ts * dt;
  //結果の表示
  if (ts % 100 == 0)
  {
    //printf("t=%f  \n",(double)t);
    if (ts >= MAX_STEP)
    { //MAX_STEPに達したらプログラムを終了
      fclose(fp);
      exit(0);
    }
  }
}

void keyboard(unsigned char key, int x, int y)
{
  double tmpcin;
  switch (key)
  {
  case 'Q':
    break;
  case 'q':
    exit(0);
    break;
  case '\033': /* '\033' = ESC */
    exit(0);
    break;
  case 'h':
    monitor.SetCenter(1.0, 0);
    break;
  case 'l':
    monitor.SetCenter(-1.0, 0);
    break;
  case 'j':
    monitor.SetCenter(0, 1.0);
    break;
  case 'k':
    monitor.SetCenter(0, -1.0);
    break;
  case 'c':
    break;
  case 'g':
    break;
  case 'z':
    monitor.SetZoom(1.1 / 1.0);
    break;
  case 'x':
    monitor.SetZoom(1.0 / 1.1);
    break;
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
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  int i, j;
  double o = 0.0, n;
  //ディスプレイ表示用文字列
  //char str[256],str1[256],str2[256],str3[256],str4[256],str5[256];
  char str[256];
  glClear(GL_COLOR_BUFFER_BIT);
  glLoadIdentity();
  //glutKeyboardFunc(keyboard);

  runge(); //runge-kuttaを回す

  monitor.SetCenter(trunk(j, 0).aa.r.get_x() / L * 2 - 1.5, 0);
  monitor.SetAllColor(0.0, 0.0, 0.0);
  for (j = 0; j < N; j++)
  {
    monitor.SetAllColor(0.0, 0.0, 0.0);
    //体節を表示
    monitor.DrawCircle(trunk(j, 0).aa.r.get_x() / L * 2 - 1.5, trunk(j, 0).aa.r.get_y() / L * 2, 0.02);
    if (j != 0 && j != N - 1)
    {

      if (F_lGround[j].get_abs() <= 0.0)
      {                                     //トルクが正であればその質点は青色
        monitor.SetAllColor(0.0, 0.0, 1.0); //青
      }
      else
      {
        monitor.SetAllColor(1.0, 0.0, .00); //赤
      }
      monitor.DrawCircle(leg(j, 0).aa.r.get_x() / L * 2 - 1.5, leg(j, 0).aa.r.get_y() / L * 2, 0.015);
    }
  }

  //リンクを描画
  monitor.SetAllColor(0.0, 0.0, 0.0);
  for (j = 0; j < N; j++)
  {
    if (j != N - 1)
    {
      //体節間のリンクを描画
      monitor.DrawLine(trunk(j, 0).aa.r.get_x() / L * 2 - 1.5, trunk(j, 0).aa.r.get_y() / L * 2, trunk(j + 1, 0).aa.r.get_x() / L * 2 - 1.5, trunk(j + 1, 0).aa.r.get_y() / L * 2, 0.005);
      if (j != 0)
      {
        //体節-脚間のリンクを描画
        monitor.DrawLine(trunk(j, 0).aa.r.get_x() / L * 2 - 1.5, trunk(j, 0).aa.r.get_y() / L * 2, leg(j, 0).aa.r.get_x() / L * 2 - 1.5, leg(j, 0).aa.r.get_y() / L * 2, 0.005);
      }
    }
  }
  //地面
  monitor.SetAllColor(0.0, 0.0, 0.0);
  monitor.DrawLine(0.3 * -50, 0.0, 0.3 * (300), 0.0, 0.005);
  for (j = -N / 3; j < 300; j++)
  {
    monitor.DrawLine(1.0 * j, 2, 1.0 * (j), -2, 0.0005);
  }
  monitor.SetAllColor(0.0, 0.0, 0.0);

  //***********力の可視化************/
  double base_x = 0.0, base_y = 0.0, F_x = 0.0, F_y = 0.0;
  Vector2D F;
  Vector2D g;
  g.set_vec(0.0, -9.8);
  //printf("GroundF_x=%lf  GroundF_y=%lf  mg_y=%lf\n",F_lGround[1].get_x(),F_lGround[1].get_y());

  ////脚
  //合力
  /*
  i = 1;
  base_x = leg(i, 0).aa.r.get_x() / L * 2 - 1.5;
  base_y = leg(i, 0).aa.r.get_y() / L * 2;
  F = (F_lBody[i] + F_lTorque[i] + F_lFriction[i] + F_lGround[i] + g * m) * (1 / m);
  F_x = base_x + F.get_x();
  F_y = base_y + F.get_y();
  monitor.SetAllColor(0.0, 0.0, 0.0);
  monitor.DrawLine(base_x, base_y, F_x, F_y, 5.0);
  //摩擦力
  monitor.SetAllColor(1.0, 0.0, 0.0);
  F = (F_lFriction[i]) * (1 / m);
  F_x = base_x + F.get_x() * 1000000000.0;
  F_y = base_y + F.get_y() * 1000000000.0;
  monitor.DrawLine(base_x, base_y, F_x, F_y, 1.5);
  //トルク
  monitor.SetAllColor(0.0, 1.0, 0.0);
  F = (F_lTorque[i]) * (1 / m);
  F_x = base_x + F.get_x();
  F_y = base_y + F.get_y();
  monitor.DrawLine(base_x, base_y, F_x, F_y, 0.5);
  //床からの反力
  monitor.SetAllColor(0.0, 0.0, 1.0);
  F = (F_lGround[i]) * (1 / m);
  F_x = base_x + F.get_x();
  F_y = base_y + F.get_y();
  monitor.DrawLine(base_x, base_y, F_x, F_y, 0.5);
  //リンクのバネダンパ
  monitor.SetAllColor(0.0, 1.0, 1.0);
  F = (F_lBody[i]) * (1 / m);
  F_x = base_x + F.get_x();
  F_y = base_y + F.get_y();
  monitor.DrawLine(base_x, base_y, F_x, F_y, 0.5);
  */
  ////体節
  //合力
  /*
  i = 1;
  base_x = trunk(1, 0).aa.r.get_x() / L * 2 - 1.5;
  base_y = trunk(1, 0).aa.r.get_y() / L * 2;
  F = (F_tBody[i] + F_tTorque[i] + F_tFriction[i] + F_tGround[i] + g * m) * (1 / m);
  F_x = base_x + F.get_x();
  F_y = base_y + F.get_y();
  monitor.SetAllColor(0.0, 0.0, 0.0);
  monitor.DrawLine(base_x, base_y, F_x, F_y, 1.0);
  //摩擦力
  monitor.SetAllColor(1.0, 0.0, 0.0);
  F = (F_tFriction[i]) * (1 / m);
  F_x = base_x + F.get_x()*1000000000.0;
  F_y = base_y + F.get_y()*1000000000.0;
  monitor.DrawLine(base_x, base_y, F_x, F_y, 5.5);
  //トルク
  monitor.SetAllColor(0.0, 1.0, 0.0);
  F = (F_tTorque[i]) * (1 / m);
  F_x = base_x + F.get_x();
  F_y = base_y + F.get_y();
  monitor.DrawLine(base_x, base_y, F_x, F_y, 0.5);
  //床からの反力
  monitor.SetAllColor(0.0, 0.0, 1.0);
  F = (F_tGround[i]) * (1 / m);
  F_x = base_x + F.get_x();
  F_y = base_y + F.get_y();
  monitor.DrawLine(base_x, base_y, F_x, F_y, 0.5);
  //リンクのバネダンパ
  monitor.SetAllColor(0.0, 1.0, 1.0);
  F = (F_tBody[i]) * (1 / m);
  F_x = base_x + F.get_x();
  F_y = base_y + F.get_y();
  monitor.DrawLine(base_x, base_y, F_x, F_y, 0.5);
*/

  monitor.SetAllColor(0.0, 0.0, 0.0);
  double x = 0.1, y = 0.85, dy = 0.1;
  sprintf(str, "Isotropic friction");
  monitor.String(x, y, str);
  sprintf(str, "ts=%d", ts);
  y = y - dy;
  monitor.String(x, y, str);
  y = y - dy;
  sprintf(str, "N=%d", N);
  monitor.String(x, y, str);
  y = y - 3.0 * dy;
  sprintf(str, "Black:TotalF");
  monitor.String(x, y, str);
  y = y - dy;
  sprintf(str, "Red:Friction");
  monitor.String(x, y, str);
  y = y - dy;
  sprintf(str, "Green:Torque");
  monitor.String(x, y, str);
  y = y - dy;
  sprintf(str, "Blue:GroundF");
  monitor.String(x, y, str);
  y = y - dy;
  sprintf(str, "Sky:SpringDamper");
  monitor.String(x, y, str);

  //力とトルクを可視化
  y = y - 5.5 * dy;
  x = -0.9;
  sprintf(str, "x0_lx=%lf,v0_lx=%lf  y0_lx=%lf,v0_ly=%lf\n", leg(1, 0).aa.r.get_x(), leg(1, 0).aa.v.get_x(), leg(1, 0).aa.r.get_y(), leg(1, 0).aa.v.get_y());
  monitor.String(x, y, str);
  y = y - dy;
  sprintf(str, "phi_t1=%lf,  phi_l1=%lf\n", phi_t[1], phi_l[1]);
  monitor.String(x, y, str);
  y = y - dy;
  sprintf(str, "Torque0_tx=%lf,  Torque0_ty=%lf\n", F_tTorque[0].get_x(), F_tTorque[0].get_y());
  monitor.String(x, y, str);
  y = y - dy;
  sprintf(str, "phi_osci[1]=%lf  sin1=%lf phi_lBar[1]=%lf,lbar_osci[1]=%lf\n", phi_osci(1, 0), sin(phi_osci(1, 0)), phi_lBar[1], lbar_osci[1]);
  monitor.String(x, y, str);

  //パラメータ値を可視化
  monitor.SetAllColor(0.0, 0.0, 0.0);
  x = 0.6, y = 0.85, dy = 0.07;
  sprintf(str, "N=%d", N);
  monitor.String(x, y, str);
  y = y - dy;
  sprintf(str, "m=%.2lf", m);
  monitor.String(x, y, str);
  y = y - dy;
  sprintf(str, "dt=%lf", dt);
  monitor.String(x, y, str);
  y = y - dy;
  sprintf(str, "omega=%.2lf", omega);
  monitor.String(x, y, str);
  y = y - dy;
  sprintf(str, "sigma=%.2lf", sigma);
  monitor.String(x, y, str);
  y = y - dy;
  sprintf(str, "Bstan=%.2lf", Bstan);
  monitor.String(x, y, str);
  y = y - dy;
  sprintf(str, "Bswin=%.2lf", Bswin);
  monitor.String(x, y, str);
  y = y - dy;
  sprintf(str, "alpha=%.2lf", alpha);
  monitor.String(x, y, str);
  y = y - dy;
  sprintf(str, "l0_osci=%.2lf", l0_osci);
  monitor.String(x, y, str);
  y = y - dy;
  sprintf(str, "kt=%.2lf", kt);
  monitor.String(x, y, str);
  y = y - dy;
  sprintf(str, "ct=%.2lf", ct);
  monitor.String(x, y, str);
  y = y - dy;
  sprintf(str, "kl=%.2lf", kl);
  monitor.String(x, y, str);
  y = y - dy;
  sprintf(str, "cl=%.2lf", cl);
  monitor.String(x, y, str);
  y = y - dy;
  sprintf(str, "kt_pas=%.2lf", kt_pas);
  monitor.String(x, y, str);
  y = y - dy;
  sprintf(str, "kt_act=%.2lf", kt_act);
  monitor.String(x, y, str);
  y = y - dy;
  sprintf(str, "kl_pas=%.2lf", kl_pas);
  monitor.String(x, y, str);
  y = y - dy;
  sprintf(str, "kl_act=%.2lf", kl_act);
  monitor.String(x, y, str);
  y = y - dy;
  sprintf(str, "ct_pas=%.2lf", ct_pas);
  monitor.String(x, y, str);
  y = y - dy;
  sprintf(str, "cl_pas=%.2lf", cl_pas);
  monitor.String(x, y, str);
      y = y - dy;
  sprintf(str, "c_gr=%.2lf",c_gr);
  monitor.String(x, y, str);
      y = y - dy;
  sprintf(str, "k_gr=%.2lf",k_gr);
  monitor.String(x, y, str);
      y = y - dy;
  sprintf(str, "mu=%.2lf",mu);
  monitor.String(x, y, str);
        y = y - dy;
  sprintf(str, "c_tanh=%.2lf",c_tanh);
  monitor.String(x, y, str);
        y = y - dy;
  sprintf(str, "l_t0=%.2lf",l_t0);
  monitor.String(x, y, str);
        y = y - dy;
  sprintf(str, "l_l0=%.2lf",l_l0);
  monitor.String(x, y, str);

  glFlush();
  glutSwapBuffers();

  //1000ステップごとに画像を保存
  if (!(ts % 100))
  {
    capture(&ts);
  }
}

void mouse(int button, int state, int x, int y)
{
  switch (button)
  {
  case GLUT_LEFT_BUTTON:
    if (state == GLUT_DOWN)
    {
      glutIdleFunc(0);
      std::cout << "left: on" << std::endl;
    }
    else
    {
      glutIdleFunc(idle);
      std::cout << "left: off" << std::endl;
    }
    break;

  case GLUT_MIDDLE_BUTTON:
    if (state == GLUT_DOWN)
    {
      glutIdleFunc(0);
      std::cout << "middle: on" << std::endl;
    }
    else
    {
      glutIdleFunc(idle);
      std::cout << "middle: off" << std::endl;
    }
    break;

  case GLUT_RIGHT_BUTTON:
    if (state == GLUT_DOWN)
    {
      glutIdleFunc(0);
      std::cout << "right: on" << std::endl;
    }
    else
    {
      glutIdleFunc(idle);
      std::cout << "right: off" << std::endl;
    }
    break;
  }
}

void resize(int w, int h)
{
  glViewport(0, 0, w, h);

  //glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  glOrtho(-w / 0.0, w / 20.0, -h / 20.0, h / 20.0, -3.0, 3.0);
  //gluPerspective( 30.0, (double)w / (double)h, 1.0 , 100.0 );
  //gluLookAt( 0.0 , 0.0 , 3.8 , 0.0 , 0.0 , 0.0 , 0.0 , 1.0 , 0.0 );
  //monitor.SetWindowSize( w , h );
  //glMatrixMode(GL_MODELVIEW);
}

void OpenGL_init(int *argcp, char **argv)
{
  init(); //初期条件を設定

  glutInit(argcp, argv);

  glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE);

  glutInitWindowSize(monitor.GetWindowSize(Monitor::X), monitor.GetWindowSize(Monitor::Y));
  glutInitWindowPosition(10, 100); //(10,100)
  winid = glutCreateWindow("simulation");
  glutDisplayFunc(display);
  glutReshapeFunc(resize);
  glutKeyboardFunc(keyboard);
  glutMouseFunc(mouse);
  glClearColor(1.0, 1.0, 1.0, 1.0);
}

void monitor_init()
{
  double zoom = 0.8;
  monitor.SetWindowSize(800, 600);
  //   monitor.SetMode( 0 );
  monitor.SetMovieMode(1);
  monitor.SetMovieName("./MovieDir/temp_");
  monitor.SetZoom(zoom);
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
  glutMainLoop(); //無限ループ

  return 0;
}
