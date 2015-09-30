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
	vector<vector<double> > quick;
	vector<vector<Nodo> >   Qnodos;
	
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
	Nodo t();
	void setQuick(int P, double a3, int E, double a1, int W1, double W1a2, int W2, double W2a2, double fu, double fv);
	vector<vector<double> > getQuick(){return quick;}
	void setQNodos(Nodo P, Nodo E, Nodo W, Nodo W1, Nodo W2);
	vector<vector<Nodo> > getQNodos(){return Qnodos;}
	
};


#endif
