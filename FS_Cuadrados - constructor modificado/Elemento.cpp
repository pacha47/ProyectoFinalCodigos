#include <iostream>

#include "Nodo.h"
#include "Arista.h"
#include "Elemento.h"

extern double dt;

/// ***************************** NS *****************************

/// ADIVINAMOS LA PRESION DEL ELEMENTO
void Elemento::setPresion( double pb){
	p = p_n = pb;}
/// LA RESTA ENTRE ELEMENTOS NO ES MAS QUE LA DIFERENCIA ENTRE SUS PRESIONES
double Elemento::operator-(Elemento e){	
	return (this->p - e.p);}

double Elemento::operator<(Elemento e){	
	return (this->p_n - e.p_n);}

/// VERIFICAMOS DE QUE LADO ESTA LA ARISTA EN EL ELEMENTO
bool Elemento::left(Nodo midA){/// si esta a la izquierda da >0
	return ((midPoint - midA)*Nodo(1,0) > .0001)? true : false; }
bool Elemento::top(Nodo midA){/// si esta en top da <0
	return ((midPoint - midA)*Nodo(0,1) < -.0001)? true : false;
}
Elemento Elemento::getVecino(int i){
	aristas[i]->getVecino(*this);
}

/// FUNCION PARA OBETENER LA ARISTA
Arista Elemento::getAri(Nodo n){
	for(int i=0;i<aristas.size();i++){
		/// HACEMOS VECTORES CON LA MISMA DIRECCION QUE LAS NORMALES
		/// Y LOS COMPARAMOS CON EL VECTOR n QUE ES UNA NORMAL ORTOGONAL
		double a = (aristas[i]->midPoint - midPoint) * n;
		if(a > 0.00001)
			return *aristas[i];
	}
	return Arista();
}

/// RETORNA LA ECUACION DE PRESION DE DICHO ELEMENTO
vec Elemento::ecuaPresion(int n, double &f){
	
	/// VARIABLES NECESARIAS
	vec ecua(n);
	int nvecino;
	/// OBTENEMOS LAS ARISTAS PARA CORREGIR LA PRESION
	Arista top = getAri(Nodo(0,1)),
		bot = getAri(Nodo(0,-1)),
		der = getAri(Nodo(1,0)),
		izq = getAri(Nodo(-1,0));
	/// DEFINIMOS EL BIJ (VELOCIDADES)
	double a=0 , aij = 0;
	f=0;
	
//	if(izq.isFront() == 2 || der.isFront() == 2 || bot.isFront() == 2 || top.isFront() == 2){
//		ecua(numero) = 1;
//		return ecua;
//	}
	
	///  **************** TOP
	Elemento veci = top.getVecino(*this);
	f-= (top.getv() + (veci.p - this->p) * dt / top.modulo ); 
	a = (top.isFront())? 0 : dt / top.modulo;
	aij += a;
	/// AGREGAMOS EL ELEMENTO VECINO A LA ECUACION
	ecua(veci.numero) = -a;
	
	///  **************** BOT
	veci =  bot.getVecino(*this);
	f+= (bot.getv() - (veci.p - this->p) * dt / bot.modulo);
	a = (bot.isFront())? 0 : dt / bot.modulo;
	aij += a;
	/// AGREGAMOS EL ELEMENTO VECINO A LA ECUACION
	ecua(veci.numero) = -a;
	
	///  **************** DER
	veci =  der.getVecino(*this);
	f-= (der.getu() + (veci.p - this->p) * dt / der.modulo );
	a = (der.isFront())? 0 : dt / der.modulo;
	aij += a;
	/// AGREGAMOS EL ELEMENTO VECINO A LA ECUACION
	ecua(veci.numero) = -a;
	
	///  **************** IZQ
	veci = izq.getVecino(*this);
	f+= (izq.getu() - (veci.p - this->p) * dt / izq.modulo);
	a = (izq.isFront())? 0 : dt / izq.modulo;
	aij += a;
	/// AGREGAMOS EL ELEMENTO VECINO A LA ECUACION
	ecua(veci.numero) = -a;
	
	/// AGREGAMOS LA VARIABLE DEL ELEMENTO CENTRAL
	ecua(this->numero) = aij;
	
	return ecua;
}
/// CORREGIMOS LA PRESION
void Elemento::corregir(double presion){
	this->p_n = this->p;
	this->p = presion;
}

/// ASIGNAMOS LAS PRESIONES A LOS NODOS DEL ELEMENTO
void Elemento::addp(){
	for(int i=0;i<nodos.size();i++){
		nodos[i]->p_aux += this->p;
		nodos[i]->cp ++;
		}
}

/// RETORNA LA ECUACION DE PRESION DE DICHO ELEMENTO
vec Elemento::setVecinos(int n){
	
	/// VARIABLES NECESARIAS
	vec ecua(n);
	/// OBTENEMOS LAS ARISTAS PARA CORREGIR LA PRESION
	Arista top = getAri(Nodo(0,1)),
		bot = getAri(Nodo(0,-1)),
		der = getAri(Nodo(1,0)),
		izq = getAri(Nodo(-1,0));
	
	ecua(top.getVecino(*this).numero) = 1.0;
	ecua(bot.getVecino(*this).numero) = 1.0;
	ecua(der.getVecino(*this).numero) = 1.0;
	ecua(izq.getVecino(*this).numero) = 1.0;
	ecua(numero) = 1.0;
	
	return ecua;
	
}

/// ************************************************************** 

/// CONSTRUCTOR
Elemento::Elemento(vector< vector<Nodo>::iterator > n, int id){
	nodos = n; numero = id;
	double nsize = n.size();
	Nodo aux = *nodos[0];
	for(int i=1;i<nsize;i++)
		aux= aux + *nodos[i];
	midPoint = aux / nsize;
	double ux = nodos[1]->x - nodos[0]->x,
		uy = nodos[1]->y - nodos[0]->y,
		vx = nodos[2]->x - nodos[0]->x,
		vy = nodos[2]->y - nodos[0]->y;
	area = (nsize == 3)? (ux*vy - vx*uy)/2.0 : (ux*vy - vx*uy);
	
	setPresion();
}

bool Elemento::operator==(Elemento e){	
	return (numero == e.numero);}


/// 			RETORNA EL ID DE LA ARISTA EN COMUN (N1,N2)
int id_ari(Nodo n1,Nodo n2){
	for(int i=0;i<n1.aris.size();i++){
		for(int j=0;j<n2.aris.size();j++){
			if(n1.aris[i] == n2.aris[j]) return n1.aris[i];
		}}
	cout<<"ERROR NO ENCONTRO ARISTA EN COMUN Y ES CUALQUI ESO"<<endl;
}


///            AGREGAMOS LAS ARISTAS
void Elemento::setAristas(vector<Arista>::iterator Ari_Begin, vector<Elemento>::iterator Ele_Begin){
	int n = nodos.size();
	for(int i=0;i<n;i++){
		aristas.push_back(Ari_Begin + id_ari(*nodos[i],*nodos[(i+1)%n]));
		aristas[i]->addElemento(Ele_Begin + numero);
	}}

/// EVALUA SI UNA ARISTA LE PERTENECE
bool Elemento::pertenece(Arista a){
	int c=0;
	for(int i=0;i<nodos.size();i++){
		if(*nodos[i] == *a.nodos[0] || *nodos[i] == *a.nodos[1] )
			c++;
	}
	return c >= 2;}

/// RETORNAMOS EL AREA
double Elemento::getArea(){
	return area;
}
