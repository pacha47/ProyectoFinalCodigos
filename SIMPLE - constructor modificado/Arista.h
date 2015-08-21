#ifndef INC_ARISTA_H
#define INC_ARISTA_H

#include <vector>
#include "vec.h"

using namespace std;

class Nodo;
class Elemento;

class Arista{
	
private:
	
public:
	
	Nodo normal, midPoint;
	vector<vector<Nodo>::iterator> nodos;
	vector<vector<Elemento>::iterator> elementos;
	
	int tipoFront, id;
	double u,v,p, modulo,
		oldV, oldU, ap;
	
	Arista(){tipoFront = 0; modulo = -1;}
	
	void addElemento(vector<Elemento>::iterator e);
	void addNodos(vector<Nodo>::iterator n1,vector<Nodo>::iterator n2);
	
	
	/// ***************************** NS *****************************
	void setFront(double u, double v);
	void setFront(double p);
	void adivinar();
	void setId(int id);
	double iterar();
	double corregir(double alpha);
	double getAp(){return ap;};
	double getv(){return v;}
	double getu(){return u;}
	void addvel();
	vec setVecinos(int n);
	vec ecuaVelocidad(int n, double &f);
	void asignarUV(double uv);
	/// ************************************************************** 
	
	bool operator==(Arista a);
	double operator*(Nodo N);
	Elemento getVecino(Elemento e);
	double getmodulo();
	int isFront(){return tipoFront;}
	Nodo getNodo(int i){return *nodos[i];}
	
};


#endif
