#ifndef INC_ARISTA_H
#define INC_ARISTA_H

#include <vector>

using namespace std;

class Nodo;
class Elemento;

class Arista{
	
private:
	
	Nodo normal, midPoint;
	vector<vector<Nodo>::iterator> nodos;
	vector<vector<Elemento>::iterator> elementos;
	
	int tipoFront, id;
	double u,v,p, modulo, dij;
	vector<double> dist;
	
public:
	/// CONSTRUSTORES Y SETEADORES
	Arista(){tipoFront = 0; modulo = -1;}
	void addElemento(vector<Elemento>::iterator e);
	void addNodos(vector<Nodo>::iterator n1,vector<Nodo>::iterator n2);
	void setFront(double u, double v);
	void setFront(double p);
	void setId(int idAri){id = idAri;}
	
	
	/// OPERADORES 
	bool operator==(Arista a);
	double operator*(Nodo N);
	double operator%(Nodo N);
	
	/// FUNCIONES BASICAS
	Elemento getVecino(Elemento e);
	Elemento getVecino(Elemento e, vector<double> &d);
	int isFront(){return tipoFront;}
	double getModulo();
	Nodo getuv();
	Nodo getn();
	double getp();
	Nodo getMidP();
	
	double d_ij(){return dij;}
	
};


#endif
