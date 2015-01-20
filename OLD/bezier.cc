#include "bezier.h"
using namespace std;

coordinate calcNormalCurve(coordinate d1, coordinate d2)
{
	coordinate n;
	n.x = g[1].x - g[0].x;
	n.y = g[1].y - g[0].y;
	n.z = g[1].z - g[0].z;
	
	if(n.z != 0)
	{
		n.z = (n.x + n.y) / n.z;
		n.x = 1;
		n.y = 1;
	}
	else
	{
		n.z = 0;
		float x = n.x;
		n.x = n.y*-1;
		n.y = x * -1;
	}
	
	return n;
}

coordinate calcCrossProd(coordinate d1, coordinate d2)
{
	coordinate n;
	
	n.x = d1.y*d2.z - d2.z*d1.y;
	n.y = d1.z*d2.x - d2.z*d1.x;
	n.z = d1.x*d2.y - d2.x*d1.y;
	
	return n;
}

coordinate interpolate(coordinate u, coordinate v, float t)
{
	coordinate r;
	r.x = u.x + t * (v.x - u.x);
	r.y = u.y + t * (v.y - u.y);
	r.z = u.z + t * (v.z - u.z);
	
	return r;
}

vector<coordinate> copyVertList(vector<vert> v)
{
	vector<coordinate> c;
	for(int i = 0; i < v.size(); i++)
	{
		coord g;
		g.x = v[i].g.x;
		g.y = v[i].g.y;
		g.z = v[i].g.z;
		
		c.push_back(g);
	}
}

coordinate calcBezier(vector<coordinate> v, int degree, float t)
{
	for(int j = 0; j <= v_degree; j++)
	{
		for(int k = 0; k < v_degree; k++)
		{
			v[k] = interpolate(v[k], v[k+1], t);
		}
	}
	
	return v[0];
}

coordinate calcuv(vector<coordinate> vert, 
				  int u, 
				  int v, 
				  float u_degree,
				  float v_degree,
				  float t)
{
	vector<coordinate> vRow; //vRow is of u_degree
	vector<coordinate> uRow; //uRow is of v_degree
	
	int uFill = 0;
	int u0 = 0;
	int v0 = 0;
	
	for(int i = 0; i <= v_degree; i++)
	{
		//clear vRow and fill it
		vRow.clear()
		for(int j = 0; j <= u_degree; j++)
		{
			vRow.push_back(vert[u0]);
			uFill++;
		}
		
		uRow.push_back(calcBezier(vRow, u_degree, t));
	}
	
	return calcBezier(uRow, v_degree, t);
}

vector<coordinate> calcPatch(vector<coordinate> vert,
							 int u_degree,
							 int v_degree,
							 int n_divisions)
{
	vector<coordinate> patch;
	float step = 1/n_divisions;
	
	for(float t = 0; t <= 1; t += step)
	{
		vector<coordinate> vertCopy = vert;
		
		for(int u = 0; u <= u_degree; u++)
		{
			for(int v = 0; v <= v_degree; v++)
			{
				coordinate b;
				b = calcuv(vertCopy, u, v, u_degree, v_degree, t);
				patch.push_back(b);
			}
		}
	}
	
	return patch;
}

vector<coordinate> calcCurve(vector<coordinate> vert,
							 int degree
							 int n_divisions)
{
	vector<coordinate> curve;
	float step = 1/n_divisions;
	
	for(float t = 0; t <= 1; t += step)
	{
		vector<coordinate> vertCopy = vert;
		for(int i = 0; i <= degree; i++)
		{
			coordinate b;
			b = calcBezier(vert, degree, t);
			curve.push_back(b);
		}
	}
	
	return curve;
}