#ifndef INC_ELEMENTO_H
#define INC_ELEMENTO_H

#include <vector>

struct vec;
class Nodo;
class Arista;

using namespace std;

class Elemento{
	
private:
	Nodo midPoint;
	double area, p,pn,
		u,v, u_12, v_12,
		u_1,v_1;
	int numero;
	vector< vector<Arista>::iterator > aristas;
	vector< vector<Nodo>::iterator > nodos;
	
public:
	
	/// 	FUNCIONES N-S FRACTIONA STEP
	void operador_CD(double dt);
	void operador_CD_AB(double dt);
	vec operador_P1(int ne, double dt, double &f);
	void setp(double preseion);
	
	double operador_P2(double dt);
	
	void setVelNodos();
	
	/// CONSTRUCTORES 
	Elemento(){}
	Elemento(vector< vector<Nodo>::iterator > n, int id);
	void setAristas(vector<Arista>::iterator Ari_Begin, vector<Elemento>::iterator Ele_Begin);

	/// OPERADORES
	bool operator==(Elemento e);
	bool operator!=(Elemento e);
	
	double getArea();
};



#endif
