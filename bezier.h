#ifndef bezier_h
#define bezier_h

#include <vector>
#include <iostream>
#include <string>
using namespace std;

/*
int render_bezier_curve(const string & vertex_type, 
						int degree, 
						const float * vertex)
{
	return 0;
}

int render_bezier_patch(const string &vertex_type,
						int u_degree, 
						int v_degree,
						const float * vertex)
{
	return 0;
}
*/

//New coordinate structure to handle everything
struct coord
{
	//geometry
	float x;
	float y;
	float z;
	
	//color
	float r;
	float g;
	float b;
	
	//normal
	float nx;
	float ny;
	float nz;
	
	//texture
	float t1;
	float t2;
	
	//set all values to 0
	void init();
};

/*
	P
	N D
	C
	T
*/					   
void setAttributeIndices(string vType, bool &g, bool &c, bool &n, bool &t);
void setAttributeIndices(string vType,
						 int &vSize,
						 int &g,
						 int &c,
						 int &n,
						 int &t);

//takes in an array of floats and converts it to a vector of coords
vector<coord> arr2vect(const float arr[], int degree, string vType);
vector<coord> arr2vect(const float arr[], 
					   int u_degree, 
					   int v_degree, 
					   string vType);

/*
	Functions that deal with coordinate computation
*/

//interpolates between two coords
coord interpolate(coord a, coord b, float t);
//multiplies a coordinate by a constant
coord mult(coord c, float m);
//subtracts two coordinates
coord subtract(coord b, coord a);

/*
	Functions that deal with bezier curve calculations
	
	The function is overloaded
		-One function will return the plotting points
		-One function will return a coordinate at time t
		 This function is used in patch calculation as well
*/

//calculates coordinates on a bezier curve given control points
vector<coord> calcBezierCurve(vector<coord> cPts, int degree, int n);
//calculates a single coordinate of a bezier curve at time t
coord calcBezierCurve(vector<coord> cPts, int degree, float t);

/*
	Functions that deal with bezier patch calculations
*/

//calculates the coordinates on a bezier patch given control points
vector<coord> calcBezierPatch(vector<coord> cPts, 
							  int u_degree, 
							  int v_degree,
							  int n);
coord calcBezierPatch(vector<coord> cPts,
					  int u_degree,
					  int v_degree,
					  float t1,
					  float t2);

/*
	Functions deal with calculating the derivative of a bezier curve/patch
*/

//generates the control points of the derivative of a bezier curve
vector<coord> getDerBezierCurve(vector<coord> cPts, int degree);
//generates the control points of the derivative of a bezier patch
void getDerBezierPatch(vector<coord> cPts,
					   vector<coord> &du,
					   vector<coord> &dv, 
					   int u_degree,
					   int v_degree,
					   int n);

#endif