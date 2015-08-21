
#include "Nodo.h"
#include <cmath>

using namespace std;

/// CONSTRUCTOR
Nodo::Nodo(double posX, double posY){
	x = posX; y = posY; modulo = sqrt(x*x + y*y); modulo = (modulo==0)? 1: modulo;
	cu=cv=cp=u=v=p=p_aux=.0;
}
Nodo::Nodo(){
	x=0;y=0;modulo=1;cu=cv=cp=u=v=p=p_aux=.0;}
/// ************************** OPERADORES ********************
/// OPERATOR ==
bool Nodo::operator==(Nodo n){
	return (x == n.x && y == n.y); }
/// RESTA ENTRE NodoS
Nodo Nodo::operator-(Nodo n){
	return Nodo(x-n.x,y-n.y);}
/// SUMA ENTRE NodoS
Nodo Nodo::operator+(Nodo n){
	return Nodo(x+n.x,y+n.y);}
/// DIVIDIMOS EN Nodo POR UN ESCALAR
Nodo Nodo::operator/(double n){
	return Nodo(x/n,y/n);}
/// PRODUCTO PUNTO ENTRE NodoS (VECTORES)
double Nodo::operator*(Nodo n){
	return (x*n.x + y*n.y);}
Nodo Nodo::operator*(double a){
	return Nodo(x*a,y*a);}
/// CALCULAMOS LA NORMAL AL SEGMENTO QUE UNE A DOS NodoS
Nodo Nodo::operator%(Nodo n){
	return Nodo( this->y - n.y, this->x - n.x);}
/// ************************** OPERADORES ********************

/// NORMALIZAMOS AL Nodo-VECTOR
void Nodo::normalize(){
	if (modulo == 1)
		return;
	double d = x*x + y*y;
	if(modulo*modulo > d+0.01 && modulo*modulo < d-0.01)
		modulo = sqrt(d);
	x = x/modulo; y = y/modulo;
	modulo = x*x + y*y;
	}

/// AGREGAMOS LA VELOCIDAD
void Nodo::addVel(double u, double v, double d){
	this->u+=(u/d); this->v+=(v/d); cv+=(1.0/d);
}
/// RETORNAMOS EL MODULO DE LA VELOCIDAD
double Nodo::modVel(){
	return sqrt(u*u + v*v);
}
void Nodo::push_ariId(int id){ aris.push_back(id);}

