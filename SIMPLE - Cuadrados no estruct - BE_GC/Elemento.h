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
	double area, p,p_prima, ap,
		u,v;
	int numero;
	vector< vector<Arista>::iterator > aristas;
	vector< vector<Nodo>::iterator > nodos;
	
public:
	
	/// 	FUNCIONES N-S FRACTIONA STEP
	vector<vec> operador_CD(int ne, double dt, double &fu, double &fv);
	vec operador_P1(int ne, double dt, double &f);
	void setp_prima(double p_Prima);
	void setuv(double u, double v);
	
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
