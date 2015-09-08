
#include "mat.h"
#include "vec.h"

/// ******************************    CONSTRUCTOR

mat::mat(int n){
	this->n = n;
	std::vector<double> r;
	for(int i=0;i<n;i++) r.push_back(.0);
	for(int i=0;i<n;i++) m.push_back(r);
}

mat::mat(int n,int m){
	this->n = n;
	std::vector<double> r;
	for(int i=0;i<m;i++) r.push_back(.0);
	for(int i=0;i<n;i++) this->m.push_back(r);
}

/// ******************************    FUNCIONES
#include<iostream>
void mat::gradConjugado(vec b, vec &x){
	
	/// SETEAMOS LAS VARIABLES DEL METODO r x v;
	vec g = b *-1.0 , p = b; /// vector direccion y residual
	double alpha, beta , pAp;
	int i=0;
	
	while(true){
		/// GUARDAMOS VALORES DE OPERACIONES QUE SE REPITEN
		vec Avk = ((*this) * p);
		
		pAp = p * Avk;  /// A * v_k
		
		/// VARIABLE alpha
		alpha = - (g*p) / pAp; 
		
		/// ACTUALIZAMOS EL VECTOR SOLUCION Y EL RESIDUO
		x = x + p*alpha;
		g = ((*this) * x) - b;
		
		if(g*g < 1e-10 ) return;
		
		if(i++ > m.size()){
			std::cout<<"NO CONVERGE GRADIENTE CONJUGADO"<<std::endl;
			return;
		}
		
		beta = (g * Avk) / pAp;
		p = p * beta - g;
	}
}

void mat::GausSeidel(vec b, vec &x){
	
	double k = .0, ap = .0;
	vec error(n);
	
	while(true){
		for(int i=0;i<n;i++){
			for(int j=0; ceros[i][j] != -1 ; j++){
				if(i != ceros[i][j]) {k-= m[i][j]*x(ceros[i][j]);}
				else                 {ap = m[i][j];}
			}
			k+=b(i);
			x(i) = k / ap;
			k=ap=0;
		}
		error = (*this) * x - b;
		if (error*error < 1e-10) return;}
}

void mat::setRow(vec r){
	std::vector<double> aux;
	int mSize = m.size(), k=0;
	for(int j=0;j<r.n;j++){
		if(r(j) != 0){
			aux.push_back(r(j));
			ceros[mSize][k++] = j;
		}
	}
	m.push_back(aux);
	ceros[mSize][k] = -1;
	this->n = mSize+1;
}

void mat::setRow(int i, vec r){
	for(int j=0;ceros[i][j] != -1 ;j++) m[i][j] = r(ceros[i][j]);
}

void mat::setSizeCeros(int n, int m){
	ceros = new int* [n];
	for(int i=0;i<n;i++)
		ceros[i] = new int[m+1];
}


/// ******************************    OPEADORES

double &mat::operator()(int i, int j){
	if(i>n-1 || j > n-1)
		return m[0][0];
	return m[i][j];
}

vec mat::operator*(vec v){
	vec m_v(n);
	for(int i=0;i<n;i++){
		for(int j=0; ceros[i][j] != -1;j++){
			m_v(i) += (m[i][j]* v(ceros[i][j]));
	}}
	return m_v;
}

mat mat::operator*(double p){
	mat m_p(n);
	for(int i=0;i<n;i++){
		for(int j=0;j<n;j++) m_p(i,j) = m[i][j] * p;}
	return m_p;
}

mat mat::operator/(double p){
	mat m_p(n);
	for(int i=0;i<n;i++){
		for(int j=0;j<n;j++) m_p(i,j) = m[i][j] / p; }
	return m_p;
}

mat mat::operator+(mat m){
	if(this->n != m.n)
		return mat(n);
	mat m_m(n);
	for(int i=0;i<n;i++){
		for(int j=0;j<n;j++){
			m_m(i,j) = (this->m[i][j] + m(i,j));
		}}
	return m_m;
}

mat mat::operator-(mat m){
	if(this->n != m.n)
		return mat(n);
	mat m_m(n);
	for(int i=0;i<n;i++){
		for(int j=0;j<n;j++){
			m_m(i,j) = (this->m[i][j] - m(i,j));
		}}
	return m_m;
}





