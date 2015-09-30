#ifndef INC_NODO_H
#define INC_NODO_H

#include<vector>
using namespace std;

class Nodo{
	
public:
	double x,y, modulo, u,v, p, c, u_aux, v_aux, p_aux;
	vector<int> aris;
	int front=0;
	
	/// CONSTRUCTORES
	Nodo(double posX, double posY);
	Nodo();
	
	/// FUNCIONES BASICAS
	void adduvp(double u, double v, double p);
	void adduvp(double u, double v, double p, int alpha);
	void defuvp();
	void push_ariId(int id);
	
	/// OPERADORES
	bool operator==(Nodo n);
	void operator=(Nodo n);
	Nodo operator-(Nodo n);
	Nodo operator+(Nodo n);
	Nodo operator/(double b);
	double operator*(Nodo n);
	Nodo operator*(double a);
	
};

#endif
