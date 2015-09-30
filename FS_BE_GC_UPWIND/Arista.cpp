#include <cmath>

#include "Nodo.h"
#include "Elemento.h"
#include "Arista.h"

#include "vec.h"
#include "mat.h"


using namespace std;

///                    SETEADORES CONSTRUCTORES
void Arista::addElemento(vector<Elemento>::iterator e){
	
	double d = sqrt( (this->midPoint - e->midP()) * (this->midPoint - e->midP()) ) ;
	elementos.push_back(e); 
	dist.push_back(d);
	
	if(tipoFront != 0){
		dij = this->modulo / ( (this->midPoint - e->midP()) * this->normal);
		if(dij < 0){
			this->normal.x *= -1; this->normal.y *= -1;
			dij*=-1;
		}
	}
	
	if(elementos.size() > 1){
		Nodo d12 = elementos[0]->midP() - elementos[1]->midP();
		d = sqrt( d12 * d12 );
		dist.push_back(d);
		dij = fabs ( this->modulo / (this->normal * d12) );
	}
}


void Arista::addNodos(vector<Nodo>::iterator n1,vector<Nodo>::iterator n2){
	/// ASIGNAMOS LOS NODOS
	nodos.push_back(n1); nodos.push_back(n2);
	/// CALCULAMOS LA NORMAL
	normal.x = (*n1).y - (*n2).y;
	normal.y = (*n2).x - (*n1).x;
	/// NORMALIZAMOS LA NORMAL Y ASIGNAMOS EL MODULO DE LA ARISTA 
	normal.modulo =  modulo = sqrt(normal.x * normal.x + normal.y * normal.y);
	normal.x = normal.x / normal.modulo;
	normal.y = normal.y / normal.modulo;
	/// RECALCULAMOS EL MODULO == 1
	normal.modulo = normal.x * normal.x + normal.y * normal.y;
	/// CALCULAMOS EL PUNTO MEDIO
	midPoint = (*nodos[0] + *nodos[1]) / 2.0;
}

///                          SETEAMOS LOS VALORES PARA LA FRONTERA
void Arista::setFront(double u, double v){
	this->u = u;  this->v = v;  tipoFront = 1;
	
	nodos[0]->front = nodos[1]->front = 1; 
	nodos[0]->u = (fabs(nodos[0]->u) > fabs(u)) ? nodos[0]->u : u; 
	nodos[0]->v = (fabs(nodos[0]->v) > fabs(v)) ? nodos[0]->v : v;
	nodos[1]->u = (fabs(nodos[1]->u) > fabs(u)) ? nodos[1]->u : u; 
	nodos[1]->v = (fabs(nodos[1]->v) > fabs(v)) ? nodos[1]->v : v;
	}
void Arista::setFront(double p){
	this->p = p; tipoFront = 2;
	nodos[0]->front = 2; nodos[0]->p = p;
	nodos[1]->front = 2; nodos[1]->p = p;
	}

///							SETEAMOS LOS VALORES DEL QUICK PARA EL ELEMENTO QUE LLAMO
void Arista::setQuick(int P, double a3, int E, double a1, int W1, double W1a2, int W2, double W2a2, double fu, double fv){
	
	/// INICIALIZAMOS EL VECTOR DEL ELEMENTO GANADOR
	vector<double> q;
	/// ELEMENTO P
	q.push_back(P); q.push_back(a3);
	/// ELEMENTO E
	q.push_back(E); q.push_back(a1);
	/// ELEMENTO INTERPOLADO W
	q.push_back(W1); q.push_back(W1a2); q.push_back(W2); q.push_back(W2a2);
	/// FUENTE EN CASO DE TENER CONTORNO
	q.push_back(fu); q.push_back(fv);
	
	/// AGREGAMOS EL VECTOR AL QUICK DE LA ARISTA
	quick.push_back(q);
}

//aristas[i]->setQNodos(this->midPoint, E.midPoint, W, N1.midPoint , N2.midPoint);
void Arista::setQNodos(Nodo P, Nodo E, Nodo W, Nodo W1, Nodo W2){
	vector<Nodo> qn;
	qn.push_back(P); qn.push_back(E); qn.push_back(W);
	qn.push_back(W1); qn.push_back(W2);
	Qnodos.push_back(qn);
}











///                          FUNCIONES BASICAS
Elemento Arista::getVecino(Elemento e, vector<double> &d){
	
	d.clear();
	if(elementos.size() == 1){
		d.push_back(.5); d.push_back(.5);
		return e;}
	if(e == (*elementos[0])){
		d.push_back(dist[1] / dist[2]);
		d.push_back(dist[0] / dist[2]);
		return (*elementos[1]);}
	d.push_back(dist[0] / dist[2]);
	d.push_back(dist[1] / dist[2]);
	return (*elementos[0]);
}

Elemento Arista::getVecino(Elemento e){
	if(elementos.size() == 1) return e;
	if(e == (*elementos[0]))  return (*elementos[1]);
	return (*elementos[0]);
}


Nodo Arista::getuv(){ return Nodo(u,v); }
double Arista::getp(){ return p; }

Nodo Arista::getn(){ return normal; }

Nodo Arista::getMidP(){ return midPoint; }

///                        RETORNAMOS EL MODULO DE LA ARISTA
double Arista::getModulo(){
	return modulo; }

///                         OPERADORES
bool Arista::operator==(Arista a){
	bool p1 = *nodos[0] == *a.nodos[0] && *nodos[1] == *a.nodos[1],
		p2 = *nodos[1] == *a.nodos[0] && *nodos[0] == *a.nodos[1];
	return p1 || p2;}

double Arista::operator*(Nodo N){
	return normal*N;}

double Arista::operator%(Nodo N){
	double p = normal*N;
	if(p<0){
		normal.x *= -1.0; normal.y *= -1.0; p *= -1.0;}
	return p;}

Nodo Arista::t(){
	return (*nodos[0] - *nodos[1]);
}


