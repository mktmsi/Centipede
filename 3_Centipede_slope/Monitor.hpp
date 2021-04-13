#ifndef _MONITOR_HPP
#define _MONITOR_HPP

#define N_SLICES 20
const double inv_pix3 = 3.0 / M_PI ;
const double pixinv_3 = M_PI / 3.0 ;
const double inv_pix2 = 2.0 / M_PI ;
const double pixinv_2 = M_PI / 2.0 ;


class Monitor {
public:

  Monitor();
  ~Monitor();

  enum AxisID{X ,Y ,N_AXIS};
  enum RgbID{R, G, B, N_RGB};
  enum WindowID{Left, Right ,N_WIN};

  /////////////////////入力関数//////////////////////
  void SetMovieMode(int takemovie);
  void SetMovieName(char *moviename);
  void SetLength( double length );
  void SetCenter( double addx , double addy );
  void SetZoom( double zoom );
  void SetWindowSize( int x , int y );
  void SetPoint( double x , double y );
  void SetColor( double r , double g , double b , int window );
  void SetAllColor( double r , double g , double b );

  /////////////////////出力関数//////////////////////
  int GetWindowSize( int xy );
  double GetPoint( int xy );
  int GetMovieMode();

  void Vertex( double x , double y );
  void Circle( double x , double y , double radius );
  void DrawCircle( double x , double y , double radius );
  void Line( double xi , double yi , double xj , double yj );
  void DrawLine( double xi , double yi , double xj , double yj, float linewidth );
  void Rectangle( double x1 , double y1 , double x2 , double y2 ); 
  void DrawRectangle( double x1 , double y1 , double x2 , double y2 );
  void String         ( double x,double y, char *string);
 
  void CenterLine();

  void SavePPMData();

private:

  int mode;
  double deltax;
  double WindowRatio;
  double inv_length;
  double length;
  double X_Center;
  double Y_Center;
  double point[N_AXIS];
  int window[N_AXIS];
  int count;
  int takemovie;
  char moviename[256];

  double *cos_mod;
  double *sin_mod;
  double color[N_RGB];

};

#endif
