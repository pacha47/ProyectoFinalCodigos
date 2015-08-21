#ifndef MAT_H
#define MAT_H

#include<vector>

struct vec;

struct mat{
	
	mat(int n);
	mat(int n,int m);
	mat(){};
	
	std::vector<std::vector<double> > m;
	int n, **ceros;;
	
	double 	&operator()(int i, int j);
	vec 	operator*(vec &v);
	mat 	operator*(double p);
	mat 	operator/(double p);
	mat 	operator+(mat &m);
	mat 	operator-(mat &m);
	
	
	void setRow(int i, vec r);
	void setRow(vec r);
	void gauss(vec &b, vec &x);
	void gradConjugado(vec &b, vec &x);
	void GausSeidel(vec &b, vec &x);
	void setSizeCeros(int n, int m);
};






#endif
