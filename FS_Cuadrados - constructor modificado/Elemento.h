#ifndef INC_ELEMENTO_H
#define INC_ELEMENTO_H

#include <vector>
#include "vec.h"

class Nodo;
class Arista;

using namespace std;

class Elemento{
	
private:
	Nodo midPoint;
	double area, p, p_n;
	int numero;
	vector< vector<Arista>::iterator > aristas;
	vector< vector<Nodo>::iterator > nodos;
	
public:
	
	/// ***************************** NS *****************************
	void addVel(double u, double v);
	vec ecuaPresion(int n, double &f);
	void setPresion( double pb = 0);
	void addp();
	bool left(Nodo midA);
	bool top(Nodo midA);
	Elemento getVecino(int i);
	Arista getAri(Nodo n);
	void corregir(double presion);
	double getp(){return p;}
	vec setVecinos(int n);
	/// ************************************************************** 
	
	Elemento(vector< vector<Nodo>::iterator > n, int id);
	void setAristas(vector<Arista>::iterator Ari_Begin, vector<Elemento>::iterator Ele_Begin);
	double operator-(Elemento e);
	double operator<(Elemento e);
	bool operator==(Elemento e);
	
	bool pertenece(Arista a);
	double getArea();
};



#endif
