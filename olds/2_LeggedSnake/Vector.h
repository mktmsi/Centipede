 #include <iostream>
#include <cmath>


using namespace std;


class Vector2D
{
private:
    double x;
    double y;
public:
    //コンストラクタ//
    Vector2D(void) : x(0.0), y(0.0){};
    Vector2D(double init_x, double init_y) : x(init_x), y(init_y){};
    //デストラクタ//
    ~Vector2D(void) {};

    void set_x(double value){ x = value;};   //Vector2Dのx座標入力
    void set_y(double value){ y = value;};   //Vector2Dのy座標入力
    void set_vec(double value_x,double value_y);//Vector2Dに数値代入
    
    void add_x(double value) { x += value;}; //Vector2Dのx座標に()内の値を足し合わせる
    void add_y(double value) { y += value;}; //Vector2Dのy座標に()内の値を足し合わせる
    
    double get_x(void) const{ return x;};     //ベクトルのx座標を出力
    double get_y(void) const{ return y;};     //ベクトルのy座標を出力
    
    double get_abs(void) const; //ベクトルの絶対値を出力
    
    Vector2D operator + (Vector2D vec) const;   //自分にvecを足した結果を出力
    Vector2D operator - (Vector2D vec) const;   //自分からvecを引いた結果を出力
    void operator += (Vector2D vec);        //自分にvecを足す
    void operator -= (Vector2D vec);        //自分からvecを引く
    double operator * (Vector2D vec) const; //内積を出力
    double operator ^ (Vector2D vec) const; //外積を出力
    Vector2D operator * (double value) const;   //ベクトルに定数を掛ける
    Vector2D operator / (double value) const;   //ベクトルを定数で割る
    Vector2D operator % (double value) const;   //定数で割ったときの剰余
    
};

void Vector2D::set_vec(double value_x,double value_y){  //Vector2Dに数値代入
    x=value_x;
    y=value_y;
}

double Vector2D::get_abs(void) const{   //ベクトルの絶対値を出力
    return(sqrt(x*x + y*y));
}

Vector2D Vector2D::operator + (Vector2D vec)const   //自分にvecを足した結果を出力
{
    Vector2D tmp;
    tmp.x = x + vec.x;
    tmp.y = y + vec.y;
    return tmp;
}

Vector2D Vector2D::operator - (Vector2D vec)const   //自分からvecを引いた結果を出力
{
    Vector2D tmp;
    tmp.x = x - vec.x;
    tmp.y = y - vec.y;
    return tmp;
}

void Vector2D::operator += (Vector2D vec)       //自分にvecを足す
{
    x = x + vec.x;
    y = y + vec.y;
}

void Vector2D::operator -= (Vector2D vec)       //自分からvecを引く
{
    x = x - vec.x;
    y = y - vec.y;
}

double Vector2D::operator * (Vector2D vec)const //内積を出力
{
    return(x * vec.x + y * vec.y);
}

double Vector2D::operator ^ (Vector2D vec)const //外積を出力
{
    return(x * vec.y - y * vec.x);
}

Vector2D Vector2D::operator * (double value) const   //ベクトルに定数を掛ける
{
    Vector2D tmp;
    tmp.x = x * value;
    tmp.y = y * value;
    return tmp;
}

Vector2D Vector2D::operator / (double value) const   //ベクトルを定数で割る
{
    Vector2D tmp;
    tmp.x = x / value;
    tmp.y = y / value;
    return tmp;
}

Vector2D Vector2D::operator % (double value) const   //剰余
{
    Vector2D tmp;
    tmp.x = fmod(x,value);
    tmp.y = fmod(y,value);
    return tmp;
}
