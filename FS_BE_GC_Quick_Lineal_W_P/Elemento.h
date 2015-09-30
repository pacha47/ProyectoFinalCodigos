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
//	vector<double> dN;
	
public:
	
	/// 	FUNCIONES N-S FRACTIONA STEP
	vector<vec> operador_CD(int ne, double dt, double &fu, double &fv);
	vec operador_P1(int ne, double dt, double &f);
	void setp(double preseion);
	void setu12(double u);
	void setv12(double v);
	
	double operador_P2(double dt);
	
	/// FUNCIONES DE VECINDAD
	vector<int> Vecino(Elemento &P);
	vector<vec> SetVecinos(int ne);
	vec 		SetVecinosP(int ne);
	void 		setVelNodos();
	void 		SetQuick();
	
	
	/// CONSTRUCTORES 
	Elemento(){}
	Elemento(vector< vector<Nodo>::iterator > n, int id);
	void setAristas(vector<Arista>::iterator Ari_Begin, vector<Elemento>::iterator Ele_Begin);

	/// OPERADORES
	bool operator==(Elemento e);
	bool operator!=(Elemento e);
	
	double getArea();
	Nodo midP(){return midPoint;}
	
	void dib();
	void dibquick(vector<Elemento> &E);
	
};



#endif
