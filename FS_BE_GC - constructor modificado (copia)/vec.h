#ifndef VEC_H
#define VEC_H

#include<vector>

struct mat;

struct vec{
	
	vec(int n);
	vec(){};
	
	std::vector<double> v;
	int n;
	
	double &operator()(int i);
	double 	operator*(vec c);
	vec		operator+(vec v);
	vec		operator-(vec v);
	vec		operator*(double v);
	vec		operator/(double v);
	vec    	operator*(mat m);
	
};


#endif
