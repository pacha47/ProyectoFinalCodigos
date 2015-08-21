#include "Nodo.h"
#include "Arista.h"
#include "Elemento.h"

#include "vec.h"
#include "mat.h"

#include <cmath>
#include <iostream>
using namespace std;

extern double flechas, Re;

///					FUNCIONES N-S FRACTIONAL STEP
void Elemento::operador_CD(double dt){
	int n = aristas.size();
	Elemento vecino;
	Nodo ui(u,v), uj, u_f, dij;
	
	double Du=0,Dv=0,Cu=0,Cv=0,Pu=0,Pv=0;
	
	for(int i = 0; i < n ; i++){
		/// OBTENEMOS EL VECINO PARA SETEAR uj
		vecino = aristas[i]->getVecino(*this);
		switch(aristas[i]->isFront()){
		case 1:
			uj = aristas[i]->getuv(); u_f = uj;
			/// CALCULAMOS EL VECTOR DIRECCION ENTRE LOS ELEMENTOS
			dij = aristas[i]->getMidP() - midPoint;
			vecino.midPoint = aristas[i]->getMidP();
			break;
		case 2:
			uj = ui; u_f = uj;
			/// CALCULAMOS EL VECTOR DIRECCION ENTRE LOS ELEMENTOS
			dij = aristas[i]->getMidP() - midPoint;
			vecino.midPoint = aristas[i]->getMidP();
			/// HACEMOS QUE LA PRESION EN LA FRONTERA SEA NULA
			vecino.p = (this->p * -1.0) + aristas[i]->getp();
			break;
		default:
			uj = Nodo(vecino.u,vecino.v);
			u_f = (uj + ui) / 2.0; 
			/// CALCULAMOS EL VECTOR DIRECCION ENTRE LOS ELEMENTOS
			dij = vecino.midPoint - midPoint;
			break;
		}
		
		/// CALCULAMOS LA PARTE DIFUSIVA
		Du+= (uj.x - this->u) * aristas[i]->getModulo() / (*aristas[i] % dij * Re);
		Dv+= (uj.y - this->v) * aristas[i]->getModulo() / (*aristas[i] % dij * Re);
		
		/// CALCULAMOS LA PARTE CONVECTIVA
		double vn = aristas[i]->getModulo() * (*aristas[i] * u_f );
		if(vn > 0){
			Cu -= vn * this->u;
			Cv -= vn * this->v;
		}else{
			Cu -= vn * uj.x;
			Cv -= vn * uj.y;
		}
		/// AGREGAMOS LA PRESION EN EL INSTANTE N
		Pu -= (vecino.p + this->p) * aristas[i]->getModulo() * aristas[i]->getn().x / 2.0;
		Pv -= (vecino.p + this->p) * aristas[i]->getModulo() * aristas[i]->getn().y / 2.0;
	}
	
	u_12 = this->u + dt*(Du + Cu + Pu) / area;
	v_12 = this->v + dt*(Dv + Cv + Pv) / area;
}



void Elemento::operador_CD_AB(double dt){
	int n = aristas.size();
	Elemento vecino;
	Nodo ui(u,v), ui_1(u_1,v_1),uj,uj_1, u_f,uf_1, dij;
	
	double Rnu=0,Rn1u=0,Rnv=0,Rn1v=0;
	for(int i = 0; i < n ; i++){
		/// OBTENEMOS EL VECINO PARA SETEAR uj
		vecino = aristas[i]->getVecino(*this);
		switch(aristas[i]->isFront()){
		case 1:
			uj = aristas[i]->getuv();    u_f = uj;
			uj_1 = aristas[i]->getuv();  uf_1 = uj;
			/// CALCULAMOS EL VECTOR DIRECCION ENTRE LOS ELEMENTOS
			dij = aristas[i]->getMidP() - midPoint;
			vecino.midPoint = aristas[i]->getMidP();
			break;
		case 2:
			uj = ui; u_f = uj;
			uj_1 = ui_1;  uf_1 = uj_1;
			/// CALCULAMOS EL VECTOR DIRECCION ENTRE LOS ELEMENTOS
			dij = aristas[i]->getMidP() - midPoint;
			vecino.midPoint = aristas[i]->getMidP();
			/// HACEMOS QUE LA PRESION EN LA FRONTERA SEA NULA
			vecino.p = this->p * -1.0 + aristas[i]->getp();;
			vecino.pn = this->pn * -1.0 + aristas[i]->getp();;
			break;
		default:
			uj = Nodo(vecino.u,vecino.v);
			uj_1 = Nodo(vecino.u_1,vecino.v_1);
			u_f = (uj + ui) / 2.0;
			uf_1 = (uj_1 + ui_1) / 2.0;
			/// CALCULAMOS EL VECTOR DIRECCION ENTRE LOS ELEMENTOS
			dij = vecino.midPoint - midPoint;
			break;
		}
		
		/// CALCULAMOS LA PARTE DIFUSIVA
		Rnu+= (uj.x - ui.x) * aristas[i]->getModulo() / (*aristas[i] % dij * Re);
		Rnv+= (uj.y - ui.y) * aristas[i]->getModulo() / (*aristas[i] % dij * Re);
		
		Rn1u+= (uj_1.x - ui_1.x) * aristas[i]->getModulo() / (*aristas[i] % dij * Re);
		Rn1v+= (uj_1.y - ui_1.y) * aristas[i]->getModulo() / (*aristas[i] % dij * Re);
		
		/// CALCULAMOS LA PARTE CONVECTIVA
		double vn = aristas[i]->getModulo() * (*aristas[i] * u_f ),
			   vn1 = aristas[i]->getModulo() * (*aristas[i] * uf_1 );
		if(vn > 0){ Rnu -= vn * ui.x; Rnv -= vn * ui.y; }
		else{ Rnu -= vn * uj.x; Rnv -= vn * uj.y; }
		
		if(vn1 > 0){ Rn1u -= vn * ui_1.x; Rn1v -= vn * ui_1.y;}
		else{ Rn1u -= vn * uj_1.x; Rn1v -= vn * uj_1.y; }
		
		
		/// AGREGAMOS LA PRESION EN EL INSTANTE N
		Rnu -= (vecino.p + this->p) * aristas[i]->getModulo() * aristas[i]->getn().x / 2.0;
		Rnv -= (vecino.p + this->p) * aristas[i]->getModulo() * aristas[i]->getn().y / 2.0;
		
		Rn1u -= (vecino.pn + this->pn) * aristas[i]->getModulo() * aristas[i]->getn().x / 2.0;
		Rn1v -= (vecino.pn + this->pn) * aristas[i]->getModulo() * aristas[i]->getn().y / 2.0;
	}
	
	u_12 = this->u + dt*(3.0*Rnu - Rn1u) / (area*2.0);
	v_12 = this->v + dt*(3.0*Rnv - Rn1v) / (area*2.0);
}



vec Elemento::operador_P1(int ne, double dt, double &f){
	int n = aristas.size();
	vec ecuacion(ne);
	f=0;
	Elemento vecino;
	Nodo vel, dij;
	for(int i = 0; i < n ; i++){
		/// OBTENEMOS EL VECINO PARA SETEAR LA VELOCIDAD EN LA FRONTERA
		vecino = aristas[i]->getVecino(*this);
		double kdij = 0;
		switch(aristas[i]->isFront()){
		case 1:
			vel = aristas[i]->getuv();
			/// CALCULAMOS EL VECTOR DIRECCION ENTRE LOS ELEMENTOS
			dij = aristas[i]->getMidP() - midPoint;
			/// CALCULAMOS LA PARTE DIFUSIVA
			kdij =  dt / (*aristas[i] % dij );
			break;
		case 2:
			vel = Nodo(u_12,v_12);
			/// CALCULAMOS EL VECTOR DIRECCION ENTRE LOS ELEMENTOS
			dij = aristas[i]->getMidP() - midPoint;
			/// CALCULAMOS LA PARTE DIFUSIVA
			kdij =  dt / (*aristas[i] % dij );
			ecuacion(numero)+= kdij;
			vecino.p = .0;
			break;
		default:
			vel = Nodo( (vecino.u_12 + u_12) /2.0 , (vecino.v_12 + v_12)/2.0);
			/// CALCULAMOS EL VECTOR DIRECCION ENTRE LOS ELEMENTOS
			dij = vecino.midPoint - midPoint;
			/// CALCULAMOS LA PARTE DIFUSIVA
			kdij =  dt / (*aristas[i] % dij );
			/// AGREGAMOS LOS VALORES CALCULADOS A LA ECUACION
			ecuacion(numero)+= kdij;
			ecuacion(vecino.numero) -= kdij;
			break;
		}
		/// AGREGAMOS LA VELOCIDAD U^(1 + 1/2)
		f-= (*aristas[i] * vel) * aristas[i]->getModulo();
		
		/// AGREGAMOS LA PRESION EN EL INSTANTE ANTERIOR
		f-= (vecino.p - this->p) * kdij;
	}
	return ecuacion;
}


void Elemento::setp(double presion){ 
	pn = p;
	p = presion; }

double Elemento::operador_P2(double dt){
	
	int n = aristas.size();
	Elemento vecino;
	Nodo dij;
	double gradP=0;
	/// SISTEMITA PARA CALCULAR LA NUEVA VELOCIDAD
	mat M(2);
	vec F(2), v(2);
	
	double u_f, v_f;
	
	for(int i = 0; i < 2 ; i++){
		/// OBTENEMOS EL VECINO PARA SETEAR LA VELOCIDAD EN LA FRONTERA
		vecino = aristas[i]->getVecino(*this);
		
		switch(aristas[i]->isFront()){
		case 1:
			/// CALCULAMOS EL VECTOR DIRECCION ENTRE LOS ELEMENTOS
			dij = aristas[i]->getMidP() - midPoint;
			gradP = 0;
			break;
		case 2:
			vecino = aristas[2]->getVecino(*this);
			if(aristas[2]->isFront()){
				/// CALCULAMOS EL VECTOR DIRECCION ENTRE LOS ELEMENTOS
				dij = aristas[2]->getMidP() - midPoint;
				gradP = 0;
			}else{
				/// CALCULAMOS EL VECTOR DIRECCION ENTRE LOS ELEMENTOS
				dij = vecino.midPoint - midPoint;
				/// CALCULAMOS EL GRADIENTE DE PRESION
				gradP = ((vecino.p - this->p) - (vecino.pn - this->pn));
				gradP *= dt /( (*aristas[2] % dij ) * aristas[2]->getModulo());
			}
			break;
		default:
			/// CALCULAMOS EL VECTOR DIRECCION ENTRE LOS ELEMENTOS
			dij = vecino.midPoint - midPoint;
			/// CALCULAMOS EL GRADIENTE DE PRESION
			gradP = ((vecino.p - this->p) - (vecino.pn - this->pn));
			gradP *= dt /( (*aristas[i] % dij ) * aristas[i]->getModulo());
			break;
		}
		
		u_f = this->u_12;
		v_f = this->v_12;
		
		/// OBTENEMOS LA NORMAL PARA DEFINIR LOS COEFICIENTES DE LA MATRIZ
		Nodo N = (aristas[i]->isFront() == 2) ? aristas[2]->getn() : aristas[i]->getn();;
		/// SETEAMOS LOS COEFICIENTE DE LA MATRIZ Y EL TERMINO DERECHO
		M(i,0) = N.x; M(i,1) = N.y;
		F(i) = N.x*u_f + N.y*v_f - gradP;
	}
	
	M.gauss(F,v);
	
	double error = (this->u - v(0))*(this->u - v(0)) + (this->v - v(1))*(this->v - v(1));
	this->u_1 = this->u; this->v_1 = this->v;
	this->u = v(0) ;  this->v = v(1) ;
	return error;
}


void Elemento::setVelNodos(){ 
	for(int i = 0 ; i < nodos.size() ; i++) nodos[i]->adduvp(u,v,p); }

///            CONSTRUCTOR
Elemento::Elemento(vector< vector<Nodo>::iterator > n, int id){
	///ASIGNAMOS EL ID Y LOS NODOS
	nodos = n; numero = id;
	
	/// ASIGNAMOS EL VALOR DE LAS VARIABLES
	u = v = u_12 = v_12	= .0;
	
	/// CALCULAMOS EL PUNTO MEDIO
	Nodo aux = *nodos[0];
	double nsize = n.size();
	for(int i=1;i<nsize;i++)
		aux= aux + *nodos[i];
	midPoint = aux / nsize;
	/// CALCULAMOS EL AREA
	double ux = nodos[1]->x - nodos[0]->x,
		uy = nodos[1]->y - nodos[0]->y,
		vx = nodos[2]->x - nodos[0]->x,
		vy = nodos[2]->y - nodos[0]->y;
	area = (ux*vy - vx*uy)/2.0;
}

/// 			RETORNA EL ID DE LA ARISTA EN COMUN (N1,N2)
int id_ari(Nodo n1,Nodo n2){
	for(int i=0;i<n1.aris.size();i++){
		for(int j=0;j<n2.aris.size();j++){
			if(n1.aris[i] == n2.aris[j]) return n1.aris[i];
		}}
	cout<<"ERROR NO ENCONTRO ARISTA EN COMUN Y ES CUALQUI ESO"<<endl;
}

///            AGREGAMOS LAS ARISTAS
void Elemento::setAristas(vector<Arista>::iterator Ari_Begin, vector<Elemento>::iterator Ele_Begin){
	int n = nodos.size();
	for(int i=0;i<n;i++){
		aristas.push_back(Ari_Begin + id_ari(*nodos[i],*nodos[(i+1)%n]));
		aristas[i]->addElemento(Ele_Begin + numero);
	}}

///            OPERADORES 
bool Elemento::operator==(Elemento e){	
	return (numero == e.numero);}

bool Elemento::operator!=(Elemento e){	
	return (numero != e.numero);}



///            RETORNAMOS EL AREA
double Elemento::getArea(){
	return area; }
