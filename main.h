#ifndef __main_h
#define	__main_h

#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include <time.h>
#include <math.h>
#include <stack>
#include <vector>
#include <iterator>
#include <algorithm>
#include <sstream>
#include <iomanip>
#include <unordered_set>

#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

using namespace std;

class Cell{
 public:
  double elev;
  int row, col;

 Cell(double e, int r, int c): elev(e), row(r), col(c) { }
};

class Point{
public:
  float x, y, z;
  Point(float i, float j, float k): x(i), y(j), z(k) {}
};

void read_file(string file_name);
void reduce_resolution();
void get_border_ocean();
double get_value(int i, int j, vector<double> *v);
void set_value(int i, int j, vector<double> *v, double value);
void FLOOD();
void FLOOD2();
void draw_grid();
void display(void);
void keypress(unsigned char key, int x, int y);
void get_minz_maxz();
void color_ground(double elev, GLfloat *shade);
void color_water(double elev, GLfloat *shade);
void hill_shade(Point p1, Point p2, Point p3, GLfloat* shade);
GLfloat ztoscreen(GLfloat z);
GLfloat ytoscreen(GLfloat y, int num_rows);
GLfloat xtoscreen(GLfloat x, int num_cols);



#endif
