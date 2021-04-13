#include <iostream>
#include <stdio.h>
//For Linux
//#include <GL/glut.h>
//For Mac
#include <GLUT/glut.h>

#include <cmath>
#include <fstream>
#include "Monitor.hpp"

Monitor:: Monitor(){

  deltax = 0.0;
  WindowRatio = 1.0;
  length = 2.0;
  inv_length = 1.0 / length;

  window[X] = 1600;
  window[Y] = 1600;

  takemovie = 0;
  count = 0;

  sprintf(moviename, "./MovieDir/");

  cos_mod = new double[N_SLICES];
  sin_mod = new double[N_SLICES];

  for( int i = 0 ; i < N_SLICES ; i++ ){
    cos_mod[i] = cos( i*2*M_PI/N_SLICES );
    sin_mod[i] = sin( i*2*M_PI/N_SLICES );
  }

}

void Monitor::String(double x,double y, char *string)
{
    glPushAttrib(GL_ENABLE_BIT);
    glDisable(GL_LIGHTING);
    glColor3d( color[R], color[G], color[B]);
    glRasterPos3f(x, y,0.0);
    while (*string) {
        glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18, *string++);
    }
    glPopAttrib();
    
}

Monitor:: ~Monitor(){
  delete [] cos_mod;
  cos_mod = NULL;

  delete [] sin_mod;
  sin_mod = NULL;
}

/////////////////////入力関数//////////////////////
void Monitor:: SetMovieMode(int takemovie){

  this->takemovie = takemovie;

}

void Monitor:: SetMovieName(char *moviename){

  sprintf(this->moviename, moviename);

}

void Monitor:: SetCenter( double addx , double addy ){

  X_Center = addx;
  Y_Center = addy;

}

void Monitor:: SetZoom( double zoom ){

  SetLength( length * zoom );

}

void Monitor:: SetWindowSize( int x , int y ){

  window[X] = x;
  window[Y] = y;

  WindowRatio = (double)window[Y] / (double)window[X];

}

void Monitor:: SetPoint( double x , double y ){

  point[X] =  2.0 * x / window[X] - 1.0;
  point[Y] = -2.0 * y / window[Y] + 1.0;

  if( mode == 2 )
    if( point[X] > 0.0 )
      deltax = 0.5;
    else
      deltax = -0.5;

  point[X] = ( point[X] - deltax ) * length / WindowRatio + X_Center;
  point[Y] = point[Y] * length + Y_Center;

  deltax = 0.0;

}


void Monitor:: SetLength( double length ){

  this->length = length;
  this->inv_length = 1.0 / length;

}

void Monitor:: SetColor( double r , double g , double b , int window ){

  color[R] = r;
  color[G] = g;
  color[B] = b;

}

void Monitor:: SetAllColor( double r , double g , double b ){

  for( int window = 0 ; window < N_WIN ; window++ ){

  color[R] = r;
  color[G] = g;
  color[B] = b;

  }

}

//////////////////////////////////////////////////


/////////////////////出力関数//////////////////////
int Monitor:: GetWindowSize( int xy ){

  return window[xy];

}

double Monitor:: GetPoint( int xy ){

  return point[xy];

}

int Monitor:: GetMovieMode(){

  return takemovie;

}

//////////////////////////////////////////////////



void Monitor:: Vertex( double x , double y ){//描画の頂点

  glVertex2d(
	     ( x - X_Center ) * WindowRatio * inv_length + deltax,
	     ( y - Y_Center ) * inv_length
	     );

}

void Monitor:: Rectangle( double x1 , double y1 , double x2 , double y2 ){//長方形

  glBegin(GL_POLYGON);

    Vertex( x1 , y1 );
    Vertex( x2 , y1 );
    Vertex( x2 , y2 );
    Vertex( x1 , y2 );

  glEnd();

}

void Monitor:: DrawRectangle( double x1 , double y1 , double x2 , double y2 ){//長方形の描画
    glColor3d( color[R], color[G], color[B] );
    Rectangle( x1 , y1 , x2 , y2 );
}

void Monitor:: Circle( double x , double y , double radius ){//円

  glBegin(GL_POLYGON);

  for(int i=0; i< N_SLICES ; i++){
    Vertex(
	   x + radius*cos_mod[i] ,
	   y + radius*sin_mod[i] 
	   );
  }
  glEnd();

}

void Monitor:: DrawCircle( double x , double y , double radius ){//円の描画
    glColor3d( color[R], color[G], color[B] );
    Circle( x , y , radius );
}

void Monitor:: Line( double xi , double yi , double xj , double yj){//線
    glBegin(GL_LINES);
    Vertex( xi , yi );
    Vertex( xj , yj );
    glEnd();

}

void Monitor:: DrawLine( double xi , double yi , double xj , double yj, float linewidth ){//線の描画

    glColor3d( color[R], color[G], color[B] );
	glLineWidth (linewidth) ;
    Line( xi , yi , xj , yj);
}

void Monitor:: CenterLine(){//中心線

    glColor3d( 0.0 , 0.0 , 0.0 );
    glBegin(GL_LINES);

    glVertex2d( 0.0 , 1.0 );
    glVertex2d( 0.0 , -1.0 );
  
    glVertex2d( 1.0 , 0.0 );
    glVertex2d( -1.0 , 0.0 );
  
    glEnd();
}

void Monitor:: SavePPMData(){//PPMの出力

  if( takemovie ){
    char filename[256];
    sprintf(filename, "%s%05d.ppm", moviename , count);

    FILE *ppmf;
    GLubyte piximage[window[Y]][window[X]][3];

    glReadPixels(0, 0, window[X], window[Y], GL_RGB, GL_UNSIGNED_BYTE, piximage);

    ppmf=fopen(filename ,"w");
    fprintf(ppmf,"P6\n");
    fprintf(ppmf,"%d %d\n",window[X],window[Y]);
    fclose(ppmf);
    ppmf=fopen(filename,"ab");
    fputc(50,ppmf);
    fputc(53,ppmf);
    fputc(53,ppmf);
    fputc(10,ppmf);

    for(int pix_j=0; pix_j<window[Y]; pix_j++){
      for(int pix_i=0; pix_i<window[X]; pix_i++){
	for(int pix_k=0; pix_k<3; pix_k++){
	  fputc(piximage[window[Y]-1-pix_j][pix_i][pix_k],ppmf);
	}
      }
    }
    fclose(ppmf);

    count++;
  }

}
