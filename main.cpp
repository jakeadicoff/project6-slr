/*
   Ethan Zhou, Jake Adicoff

   Code to visualize flooding due to sea level rise. Given a terain,
   finds ocean after sea level has increased. Renders results in
   OpenGL.
 */

#include "main.h"
#include <chrono>

//whenever the user rotates and translates the scene, we update these
//global translation and rotation
GLfloat pos[3] = {0,0,0};
GLfloat theta[3] = {0,0,0};

// draw polygons line or filled. This will be used when rendering the
// surface.
GLint fillmode = 0;

#define ONE_COLOR 0
#define CODE_COLOR 1
#define MYCODE_COLOR 2
#define NB_COLORMAP_CHOICES 3
int COLORMAP = ONE_COLOR;

const int WINDOWSIZE = 500;
const int BIGINT = 0x0fffffff;

double max_slr;
double increment;
// for rendering
int display_slr = 0;
int density = 1;

int NODATA = -9999;
int num_rows;
int num_cols;
double minz, maxz, zrange;
stack<Cell> a;
stack<Cell> b;
vector<double> elevation_grid;
vector<double> flood_level;
unordered_set<int> visited_indices;

//for second alg
vector<stack< Cell > > flood_level_stacks;


int main(int argc, char** argv) {
  if(argc < 4 || argc > 6) {
    cout << "incorrect number of arguments:\n" <<
      "arg 1: input file name\n" <<
      "arg 2: maximum sea level rise\n" <<
      "arg 3: flooding increment\n" <<
      "arg 4: (optional) grid resolution reduction factor\n" <<
      "arg 5: (optional) flooding algorithm FLOOD or FLOOD2\n";
    exit(0);
  }

  //get user inputs
  string input_file_name = argv[1];
  max_slr = atof(argv[2]);
  increment = atof(argv[3]);
  if(argc >= 5) density = atoi(argv[4]);
  string alg;
  if(argc >= 6) alg = argv[5];

  //calculate flood levels
  read_file(input_file_name);

  auto start = std::chrono::steady_clock::now();
  if(density > 1) reduce_resolution();

  vector<double> dummy1(elevation_grid.size(),NODATA);
  flood_level = dummy1;
  get_border_ocean();

  //2 options for flooding
  if(alg == "FLOOD" || alg == "") FLOOD();
  else if(alg == "FLOOD2") {
    stack<Cell> s;
    vector<stack< Cell > > dummy2(max_slr+1,s);
    flood_level_stacks = dummy2;
    FLOOD2();
  }

  //timing
  auto end = std::chrono::steady_clock::now();
  std::chrono::duration<double> diff = end-start;
  std::cout << "Time to execute: " << diff.count() << " s\n";

  //rendering
  get_minz_maxz();
  zrange = maxz - minz;
  /* OPEN GL STUFF */
  /* open a window and initialize GLUT stuff */
  glutInit(&argc, argv);
  glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB | GLUT_DEPTH);
  glutInitWindowSize(WINDOWSIZE, WINDOWSIZE);
  glutInitWindowPosition(100,100);
  glutCreateWindow(argv[0]);

  /* register callback functions */
  glutDisplayFunc(display);
  glutKeyboardFunc(keypress);

  /* OpenGL init */
  /* set background color black*/
  glClearColor(0, 0, 0, 0);
  glEnable(GL_DEPTH_TEST);
  /* start the event handler */
  glutMainLoop();

  return 0;
}

//Populates flood_levels array, such that each index has what level it
//floods at
void FLOOD() {
  stack<Cell> *current_level = &a;
  stack<Cell> *next_level = &b;
  double current_slr = 0;

  //increment the current sea level rise value as the user requests
  for(current_slr = 0; current_slr <= max_slr; current_slr+=increment) {

    unordered_set<int> next_level_indecies;

    //run DFS to blob search all currently underwater points
    while(!current_level->empty()) {
      Cell c = current_level->top();
      current_level->pop();

      //if this point is already classified, continue
      if(get_value(c.row,c.col,&flood_level) != NODATA) continue;

      // if a cell is on the stack and has elevation that would be
      // flooded by the current sea level rise, classify it as such
      if(c.elev == NODATA || c.elev <= current_slr) {
	set_value(c.row, c.col, &flood_level, current_slr);
	for(int i = -1; i < 2; i++) {
	  for(int j = -1; j < 2; j++) {
	    // if no change or out of bounds, continue
	    if(i==0 && j==0) continue;
	    if(c.row+i > num_rows-1 || c.row + i < 0) continue;
	    if(c.col+j > num_cols-1 || c.col + j < 0) continue;

	    Cell n(get_value(c.row+i,c.col+j,&elevation_grid),c.row+i,c.col+j);

	    // if neighbor is an already classified cell, continue
	    if(get_value(n.row,n.col,&flood_level) != NODATA) continue;

	    // if neighbor is also flooded , add to current stack.
	    if(n.elev == NODATA || n.elev <= current_slr) current_level->push(n);
	    // Otherwise, add to the stack of the next level
	    // but... don't add n into the the next_level stack if
	    // it's already in there
	    else{
	      int idx = n.row*num_cols + n.col;

	      unordered_set<int>::const_iterator got =
	      	next_level_indecies.find(idx);

	      //if the index is not in the set, add it to next_level and the set
	      if (got == next_level_indecies.end()){
	      	next_level_indecies.insert(idx);
		next_level->push(n);
	      }
	    }//if/else point is above/below elevation
	  }//for j
	}//for i
      }// if(point is underwater)

      //Else, if the current point is above the sea level, add the
      //point to the next level stack
      else next_level->push(c);
    } // DFS while loop

    //swap the "current" and "next" level pointers when the current
    //stack is empty
    if(current_level == &a) {
      current_level = &b;
      next_level = &a;
    }
    else {
      current_level = &a;
      next_level = &b;
    }
  } // for current_slr
} // FLOOD


// a modification of our original alg. Runs in O(n) because a point is
// only ever pushed to any stack once and poped from a stack once.
// IE every point is pushed to the correct stack instantly
void FLOOD2() {
  flood_level_stacks[0] = a;

  for(int current_slr = 0; current_slr < max_slr; ++current_slr) {

    while(!flood_level_stacks[current_slr].empty()) {
      // pop from top, classify sea level, and determine what to do with neighbors
      Cell c = flood_level_stacks[current_slr].top();
      flood_level_stacks[current_slr].pop();

      //if the point is in a stack, it must be flooded
      set_value(c.row,c.col,&flood_level,current_slr);

      //find neighbors:
      for(int i = -1; i < 2; i++) {
	for(int j = -1; j < 2; j++) {
	  //skip degenerate cases
	  if(i==0 && j==0) continue;
	  if(c.row+i > num_rows-1 || c.row + i < 0) continue;
	  if(c.col+j > num_cols-1 || c.col + j < 0) continue;
	  //make neighbor object
	  Cell n(get_value(c.row+i,c.col+j,&elevation_grid),c.row+i,c.col+j);

	  // get index (in row major)
	  int idx = n.row*num_cols + n.col;
	  // check for existance of index
	  unordered_set<int>::const_iterator got =
	    visited_indices.find(idx);

	  //if the index is not in the set, add it to appropriate stack
	  if (got == visited_indices.end()){
	    visited_indices.insert(idx);
	    // if not flooded by current slr
	    if(n.elev > current_slr) {
	      if(n.elev <= max_slr) {
		// add to flood level that is equal to the elevation
		flood_level_stacks[n.elev].push(n);
	      } else {
		// if never flooded, add to trash stack. delete this maybe
		flood_level_stacks[max_slr].push(n);
	      }
	    } else {
	      // if found in the round of bfs, belongs in this flood level. push it
	      flood_level_stacks[current_slr].push(n);
	    }
	  }//if not visited
	}// for i
      }// for j
    }//while
  }//for current_slr
}//end func



//get all grid points that are on the border and are NODATA or 0. Push
//them onto our initial stack where we begin our DFS
void get_border_ocean() {
  for(int i = 0; i < num_cols; i++) {
    //top edge
    double top_val = get_value(0,i,&elevation_grid);
    if(top_val == NODATA || top_val >= 0) {
      Cell temp(get_value(0,i,&elevation_grid), 0, i);
      a.push(temp);
      int idx = temp.row*num_cols + temp.col;
      visited_indices.insert(idx);
    }
    //bottom edge
    double bottom_val = get_value(num_rows-1,i,&elevation_grid);
    if(bottom_val == NODATA || bottom_val >= 0) {
      Cell temp(get_value(num_rows-1,i,&elevation_grid), num_rows-1, i);
      a.push(temp);
      int idx = temp.row*num_cols + temp.col;
      visited_indices.insert(idx);
    }
  }

  for(int i = 0; i < num_rows; i++) {
    //left edge
    double left_val = get_value(i,0,&elevation_grid);
    if(left_val == NODATA || left_val >= 0) {
      Cell temp(get_value(i, 0,&elevation_grid), i, 0);
      a.push(temp);
      int idx = temp.row*num_cols + temp.col;
      visited_indices.insert(idx);
    }

    //right edge
    double right_val = get_value(i,num_cols-1,&elevation_grid);
    if(right_val == NODATA || right_val >= 0) {
      Cell temp(get_value(i, num_cols-1,&elevation_grid), i, num_cols-1);
      a.push(temp);
      int idx = temp.row*num_cols + temp.col;
      visited_indices.insert(idx);
    }
  }
}

//convenience functions
double get_value(int i, int j, vector<double> *v) {
  return v->at(i*num_cols + j);
}

void set_value(int i, int j, vector<double> *v, double value) {
  v->at(i*num_cols + j) = value;
}


//read in .asc file data
void read_file(string file_name) {
  //open file
  ifstream file_stream;
  file_stream.open(file_name.c_str(), ios::in);

  if(!file_stream.good())
    cout << "Error: not able to open file" << endl;

  string line;

  //parse file for nrows, ncols, nodata, and data
  int i = 0;
  while(getline(file_stream,line)) {
    //put grid data into row major grid
    if(i > 5) {
      //stackoverflow.com/questions/236129/most-elegant-way-to-split-a-string
      istringstream iss(line);
      vector<string> tokens{istream_iterator<string>{iss},
	  istream_iterator<string>{}};

      //stackoverflow.com/questions/20257582/convert-vectorstdstring-to-vectordouble
      vector<double> doubleVector(tokens.size());
      transform(tokens.begin(), tokens.end(), doubleVector.begin(),
		[](const std::string& val) { return std::stod(val); });

      //Append new line that was just read into file
      elevation_grid.insert(elevation_grid.end(),
			    doubleVector.begin(), doubleVector.end());
    }
    //first line is num_cols
    else if(i == 0) sscanf(line.c_str(), "%*s %d", &num_cols);
    //second is num_rows
    else if(i == 1) sscanf(line.c_str(), "%*s %d", &num_rows);
    //6th is nodata
    else if(i == 5) {
      sscanf(line.c_str(), "%*s %d", &NODATA);
      if(NODATA == 0) NODATA = -9999; //NODATA can't be 0
    }
    i++; //increment line counter
  }
}

//averages out density^2 cells into each cell
void reduce_resolution() {
  int reduced_num_rows = num_rows/density;
  int reduced_num_cols = num_cols/density;
  cout << "num rows " << num_rows << endl;
  cout << "num cols " << num_cols << endl;
  cout << "reduced rows " << reduced_num_rows << endl;
  cout << "reduced cols " << reduced_num_cols << endl;
  vector<double> elevation_grid_copy = elevation_grid;
  elevation_grid.clear();
  elevation_grid.resize(reduced_num_cols * reduced_num_rows);

  //loop over reduced array
  for(int i = 0; i < reduced_num_rows; ++i) {
    for(int j = 0 ; j < reduced_num_cols; ++j) {
      double sum = 0;
      bool is_sea_level = false;

      //loop over larger array, average into reduced array
      for(int dx = 0; dx < density; ++dx) {
	for(int dy = 0; dy < density; ++dy) {
	  double val = get_value(i*density+dx, j*density+dy, &elevation_grid_copy);
	  if(val == 0){
	    is_sea_level = true;
	    break;
	  }
	  sum += val;
	}//dy
      }//dx

      //if ANY of the squares containe a sea level cell, set the whole
      //cell to sea level
      if(is_sea_level) sum = 0;
      //note: because of integer division, we will never go out of
      //bounds on the large elevation_grid. We throw out some values
      //near the edges
      elevation_grid[i*reduced_num_cols + j] = sum/(density*density);
    }//j
  }//i

  num_rows = reduced_num_rows;
  num_cols = reduced_num_cols;
}

void keypress(unsigned char key, int x, int y) {
  switch(key) {
  case '+':
    if(display_slr + increment <= max_slr){
      display_slr += increment;
    glutPostRedisplay();
    }
    break;
  case '-':
    if(display_slr - increment >= 0) {
      display_slr -= increment;
      glutPostRedisplay();
    }
    break;
  case 'x':
    theta[0] -= 5.0;
    glutPostRedisplay();
    break;
  case 'y':
    theta[1] += 5.0;
    glutPostRedisplay();
    break;
  case 'z':
    theta[2] += 5.0;
    glutPostRedisplay();
    break;
  case 'X':
    theta[0] += 5.0;
    glutPostRedisplay();
    break;
  case 'Y':
    theta[1] -= 5.0;
    glutPostRedisplay();
    break;
  case 'Z':
    theta[2] -= 5.0;
    glutPostRedisplay();
    break;
  case '2':
    theta[0] = theta[1] = theta[2] = 0;
    glutPostRedisplay();
    break;
  case 'q':
    exit(0);
    break;
  }
}

/* this function is called whenever the window needs to be rendered */
void display(void) {
  //clear the screen
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  //clear all modeling transformations
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();

  /* The default GL window is x=[-1,1], y= [-1,1] with the origin in
     the center.  The view frustrum was set up from z=-1 to z=-10. The
     camera is at (0,0,0) looking along negative z axis.
  */

  /* First we translate and rotate our local reference system with the
    user transformation. pos[] represents the cumulative translation
    entered by the user, and theta[] the cumulative rotation entered
    by the user */
  glTranslatef(pos[0], pos[1], pos[2]);
  glRotatef(theta[0], 1,0,0); //rotate theta[0] around x-axis, etc
  glRotatef(theta[1], 0,1,0);
  glRotatef(theta[2], 0,0,1);

  /* We translated the local reference system where we want it to be; now we draw the
     object in the local reference system.  */

  draw_grid();

  //don't need to draw a cube but I found it nice for perspective
  //  cube(1); //draw a cube of size 1

  glFlush();
}

void draw_grid() {
  //draw two triangles for each grid cell. shade with hill_shade
  //dot product calculation.

  glBegin(GL_TRIANGLES);
  for (int i=0; i < num_rows-1; i += 1) {
    for (int j=0; j < num_cols-1; j += 1) {

      GLfloat shade[3];

      GLfloat p1_elev = get_value(i,j,&elevation_grid);
      GLfloat p2_elev = get_value(i+1,j,&elevation_grid);
      GLfloat p3_elev = get_value(i,j+1,&elevation_grid);
      GLfloat p4_elev = get_value(i+1,j+1,&elevation_grid);

      GLfloat p1_flood = get_value(i,j,&flood_level);
      GLfloat p2_flood = get_value(i+1,j,&flood_level);
      GLfloat p3_flood = get_value(i,j+1,&flood_level);
      GLfloat p4_flood = get_value(i+1,j+1,&flood_level);

      bool is_water = false;

      if(p1_flood <= display_slr && p2_flood <= display_slr &&
	 p3_flood <= display_slr && p4_flood <= display_slr &&
	 p1_flood != NODATA && p2_flood != NODATA &&
	 p3_flood != NODATA && p4_flood != NODATA) {
	is_water = true;
	color_water(p1_elev, shade);
	p1_elev = p2_elev = p3_elev = p4_elev = display_slr;
      } else {
	color_ground(p1_elev,shade);
      }

      //triangle 1
      //add hill shade if ground
      if(!is_water){
	Point p1(float(i+1), float(j), p2_elev);
	Point p2(float(i), float(j), p1_elev);
	Point p3(float(i), float(j+1), p3_elev);
	hill_shade(p1, p2, p3, shade);
      }

      glColor3fv(shade);
      glVertex3f(xtoscreen(i,num_cols),
		 ytoscreen(j,num_rows),
		 ztoscreen(p1_elev));

      glVertex3f(xtoscreen(i+1,num_cols),
		 ytoscreen(j,num_rows),
		 ztoscreen(p2_elev));

      glVertex3f(xtoscreen(i,num_cols),
		 ytoscreen(j+1,num_rows),
		 ztoscreen(p4_elev));

      //triangle 2
      if(!is_water){
	Point pa(float(i+1), float(j+1), p4_elev);
	Point pb(float(i+1), float(j), p1_elev);
	Point pc(float(i), float(j+1), p2_elev);

	hill_shade(pa, pb, pc, shade);
	glColor3fv(shade);
      }
      glVertex3f(xtoscreen(i,num_cols),
		 ytoscreen(j+1,num_rows),
		 ztoscreen(p4_elev));

      glVertex3f(xtoscreen(i+1,num_cols),
		 ytoscreen(j+1,num_rows),
		 ztoscreen(p3_elev));

      glVertex3f(xtoscreen(i+1,num_cols),
		 ytoscreen(j,num_rows),
		 ztoscreen(p2_elev));
    }
  }
  glEnd();
}

void color_ground(double elev, GLfloat *shade) {
  double pct_height = elev/zrange;

  float high[3] = {.52, 1, .4};
  float low[3] = {.31, .35, .3};

  shade[0] = (1 - pct_height) * low[0] + pct_height * high[0];
  shade[1] = (1 - pct_height) * low[1] + pct_height * high[1];
  shade[2] = (1 - pct_height) * low[2] + pct_height * high[2];
}

void color_water(double elev, GLfloat *shade) {
  double pct_height = elev/display_slr;

  float high[3] = {.11, .71, .88};
  float low[3] = {0, 0, .8};

  shade[0] = (1 - pct_height) * low[0] + pct_height * high[0];
  shade[1] = (1 - pct_height) * low[1] + pct_height * high[1];
  shade[2] = (1 - pct_height) * low[2] + pct_height * high[2];
}

//calculates how bright a triangle should be, based on how much it
//faces the sun. Uses the dot product between the incident sun vector
//and the normal vector of the triangle.
void hill_shade(Point p1, Point p2, Point p3, GLfloat* shade){
  //calculate normal vector from triangle
  Point U(p2.x - p1.x, p2.y - p1.y, p2.z - p1.z);
  Point V(p3.x - p1.x, p3.y - p1.y, p3.z - p1.z);
  Point N(U.y*V.z - U.z*V.y, U.z*V.x-U.x*V.z, U.x*V.y - U.y*V.x);

  //normalize N:
  float n_len = sqrt(N.x*N.x + N.y*N.y + N.z*N.z);
  N.x = N.x/n_len;
  N.y = N.y/n_len;
  N.z = N.z/n_len;

  //calculate dot product between sun vector and normal vector
  Point sun_incidence(-0.577, -0.577, -0.577);
  float dot_product =
    N.x*sun_incidence.x + N.y*sun_incidence.y + N.z*sun_incidence.z;

  float intensity = .5;
  shade[0] = shade[0]*(intensity*dot_product) + (1-intensity)*shade[0];
  shade[1] = shade[1]*(intensity*dot_product) + (1-intensity)*shade[1];
  shade[2] = shade[2]*(intensity*dot_product) + (1-intensity)*shade[2];
}


/* x is a value in [minx, maxx]; it is mapped to [-1,1] */
GLfloat xtoscreen(GLfloat x, int num_cols) {
  //return (-1 + 2*x/WINDOWSIZE);
  return (-1 + 2*(x)/float(num_cols));
}


/* y is a value in [miny, maxy]; it is mapped to [-1,1] */
GLfloat ytoscreen(GLfloat y, int num_rows) {
  return (-1 + 2*(y)/float(num_rows));
}

/* z is a value in [minz, maxz]; it is mapped so that [0, maxz] map to [0,1] */
GLfloat ztoscreen(GLfloat z) {
  return (-1 + 2*(z-minz)/(maxz-minz))/1.5;
}

void get_minz_maxz() {
  minz = maxz = 0;
  for(int i = 0; i < num_rows; i++) {
    for(int j = 0; j < num_cols; j++) {
      double e = get_value(i,j,&elevation_grid);
      if(e > maxz) {
	maxz = e;
      } else if(e < minz) {
	minz = e;
      }
    }
  }
}
