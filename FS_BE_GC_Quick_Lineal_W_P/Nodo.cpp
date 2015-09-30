
#include "Nodo.h"
#include <cmath>

#include <iostream>

using namespace std;

/// CONSTRUCTOR
Nodo::Nodo(double posX, double posY){
	x = posX; y = posY; modulo = sqrt(x*x + y*y); modulo = (modulo==0)? 1: modulo;
	c= u = v = p = u_aux = v_aux = p_aux = 0;
}
Nodo::Nodo(){
	x=0;y=0;modulo=1;u=v=c=p=0;}

/// FUNCIONES BASICAS
void Nodo::adduvp(double u, double v, double p){ 
	if(front==1){
		c++;
		this->p_aux += p;
		return;
	}
	if(front==2){
		this->u_aux += u; this->v_aux += v; c++;
		return;
	}
	this->u_aux += u; this->v_aux += v; c++; this->p_aux += p; 
}


void Nodo::defuvp(){
	if(front==1){
		p = p_aux / c;
		c = p_aux = 0;
		return;
	}
	if(front==2){
		u = u_aux / c; v = v_aux / c; 
		c = u_aux = v_aux = 0;
		return;
	}
	u = u_aux / c; v = v_aux / c; p = p_aux / c; c = u_aux = v_aux = p_aux = 0;}

void Nodo::push_ariId(int id){ aris.push_back(id);}

/// ************************** OPERADORES ********************
void Nodo::operator=(Nodo n){
	x = n.x; y = n.y; modulo = sqrt(x*x + y*y);}
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
/// ************************ FIN - OPERADORES ******************

