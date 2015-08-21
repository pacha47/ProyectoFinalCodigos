#ifndef INC_NODO_H
#define INC_NODO_H

#include<vector>
class Nodo{
	
public:
	double x,y, modulo, u,v,p ,p_aux,cu,cv,cp;
	std::vector<int> aris;
	
	Nodo(double posX, double posY);
	Nodo();
	bool operator==(Nodo n);
	Nodo operator-(Nodo n);
	Nodo operator+(Nodo n);
	Nodo operator/(double b);
	double operator*(Nodo n);
	Nodo operator*(double a);
	Nodo operator%(Nodo n); /// retorna la normal al segmento que une los Nodos
	
	void addVel(double u, double v, double d);
	void getVel(double &u, double &v);
	double norvel(){u/=cu; v/=cv; cu=cv=0; return modVel();}
	double norp(){p = p_aux / cp; p_aux = cp = 0; return p;}
	void push_ariId(int id);
	
	double modVel();
	void normalize();
	
};

#endif
