#ifndef BEZ_FUNC_H
#define BEZ_FUNC_H

#include <vector>
using namespace std;

struct coordinate
{
	float x;
	float y;
	float z;
};

struct vert
{
	coordinate g;
	coordinate c;
	coordinate n;
	coordinate t;
};

//calculate normal vector based on two coordinates
coordinate calcNormalCurve(coordinate d1, coordinate d2);
coordinate calcCrossProd(coordinate d1, coordinate d2);
//interpolates between two points at some time t
coordinate interpolate(coordinate u, coordinate v, float t);
//copies a vert list and outputs a coordinate list
vector<coordinate> copyVertList(vector<vert> v);
//calculates a bezier curve at time t given control points in t
coordinate calcBezier(vector<coordinate> v, int degree, float t);
//calculates a bezier patch at time t given control points v
coordinate calcuv(vector<coordinate> vert,
				  int u,
				  int v,
				  int u_degree,
				  int v_degree,
				  float t);
//outputs a vector of coordinates of the bezier patch
vector<coordinate> calcPatch(vector<coordinate> vert,
							 int u_degree,
							 int v_degree
							 int n_divisions);
//outputs a vector of coordinates of the bezier curve
vector<coordinate> calcCurve(vector<coordinate> vert,
							 int degree
							 int n_divisions);

#endif 