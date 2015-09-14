
#include "vec.h"
#include "mat.h"

/// ******************************    CONSTRUCTOR
vec::vec(int n){
	this->n = n;
	for(int i=0;i<n;i++) v.push_back(.0);
}

/// ******************************    OPEADORES
double &vec::operator()(int i){
	if(i>n-1)
		return v[0];
	return v[i];
}

double vec::operator*(vec c){
	if(this->n != c.n)
		return 0;
	double s=0;
	for(int i=0;i<c.n;i++) s+= (this->v[i] * c(i));
	return s;
}


vec vec::operator+(vec c){
	if(this->n != c.n)
		return vec(n);
	vec v_sum(n);
	for(int i=0;i<c.n;i++) v_sum(i) = v[i] + c(i);
	return v_sum;
}

vec vec::operator-(vec c){
	if(this->n != c.n)
		return vec(n);
	vec v_res(n);
	for(int i=0;i<c.n;i++) v_res(i) = v[i] - c(i);
	return v_res;
}


vec vec::operator/(double c){
	vec v_res(n);
	for(int i=0;i<n;i++) v_res(i) = v[i] / c;
	return v_res;
}


vec vec::operator*(double c){
	vec v_res(n);
	for(int i=0;i<n;i++) v_res(i) = v[i] * c;
	return v_res;
}

vec vec::operator*(mat m){
	if(m.n != n) return vec(n);
	vec v_m(n);
	for(int i=0;i<n;i++){
		for(int j=0;j<n;j++){
			v_m(i) += v[j] * m(j,i);
		}
	}
	return v_m;
}
