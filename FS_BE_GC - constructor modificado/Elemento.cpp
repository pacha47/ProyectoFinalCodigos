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
vector<vec> Elemento::operador_CD(int ne, double dt, double &fu, double &fv){
	int n = aristas.size();
	Elemento vecino;
	Nodo ui(u,v), uj, u_f, dij;
	
	vec ecuau(2*ne);
	vec ecuav(2*ne);
	double aP=0,D=0,C=0,Pu=0,Pv=0;
	fu = fv = 0;
	vector<double> dist;
	
	for(int i = 0; i < n ; i++){
		/// OBTENEMOS EL VECINO PARA SETEAR uj
		vecino = aristas[i]->getVecino(*this, dist);
		switch(aristas[i]->isFront()){
		case 1:
			uj = aristas[i]->getuv(); u_f = uj;
			
			fu += uj.x * aristas[i]->d_ij() / Re;
			fv += uj.y * aristas[i]->d_ij() / Re;
			
			break;
		case 2:
			uj = ui; u_f = uj;
			/// HACEMOS QUE LA PRESION EN LA FRONTERA SEA NULA
			vecino.p = (this->p * -1.0) + aristas[i]->getp();
			break;
		default:
			uj = Nodo(vecino.u,vecino.v);
			u_f = uj * dist[1] + ui * dist[0]; 
			/// CALCULAMOS EL VECTOR DIRECCION ENTRE LOS ELEMENTOS
			*aristas[i] % (vecino.midPoint - midPoint);
			break;
		}
		
		/// CALCULAMOS LA PARTE DIFUSIVA
		D = aristas[i]->d_ij() / Re;
		aP+=D;
		
		/// CALCULAMOS LA PARTE CONVECTIVA
		double vn = aristas[i]->getModulo() * (*aristas[i] * u_f );
		if(vn > 0){
			aP += vn;
			C = 0;
		}else{
			C = vn;
		}
		/// AGREGAMOS LA PRESION EN EL INSTANTE N
		Pu -= (vecino.p * dist[1] + this->p * dist[0] ) * aristas[i]->getModulo() * aristas[i]->getn().x;
		Pv -= (vecino.p * dist[1] + this->p * dist[0] ) * aristas[i]->getModulo() * aristas[i]->getn().y;
		
		ecuau(vecino.numero*2  ) = C - D;
		ecuav(vecino.numero*2+1) = C - D;
	}
	
	ecuau(this->numero*2  ) = aP + this->area / dt;
	ecuav(this->numero*2+1) = aP + this->area / dt;
	
	fu += this->u * this->area / dt + Pu;
	fv += this->v * this->area / dt + Pv;
	
	vector<vec> ecuacion;
	ecuacion.push_back(ecuau); ecuacion.push_back(ecuav); 
	return ecuacion;
}


vec Elemento::operador_P1(int ne, double dt, double &f){
	int n = aristas.size();
	vec ecuacion(ne);
	f=0;
	Elemento vecino;
	Nodo vel, dij;
	vector<double> dist;
	
	for(int i = 0; i < n ; i++){
		/// OBTENEMOS EL VECINO PARA SETEAR LA VELOCIDAD EN LA FRONTERA
		vecino = aristas[i]->getVecino(*this, dist);
		double kdij = dt * aristas[i]->d_ij() / aristas[i]->getModulo();
		switch(aristas[i]->isFront()){
		case 1:
			vel = aristas[i]->getuv();
			break;
		case 2:
			vel = Nodo(u_12,v_12);
			ecuacion(numero)+= kdij;
			vecino.p = .0;
			break;
		default:
			vel = Nodo(vecino.u_12 * dist[1] + u_12 * dist[0], vecino.v_12 * dist[1] + v_12 * dist[0]);
			/// AGREGAMOS LOS VALORES CALCULADOS A LA ECUACION
			ecuacion(numero)+= kdij;
			ecuacion(vecino.numero) -= kdij;
			break;
		}
		*aristas[i] % (vecino.midPoint - midPoint);
		/// AGREGAMOS LA VELOCIDAD U^(1 + 1/2)
		f-= (*aristas[i] * vel) * aristas[i]->getModulo();
		
		/// AGREGAMOS LA PRESION EN EL INSTANTE ANTERIOR
		f-= (vecino.p - this->p) * kdij;
	}
	return ecuacion;
}

void Elemento::setu12(double u){ 
	this->u_12 = u;
}

void Elemento::setv12(double v){ 
	this->v_12 = v;
}

void Elemento::setp(double presion){ 
	pn = p;
	p = presion; }

double Elemento::operador_P2(double dt){
	
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
	mat M(2);
	vec F(2), v(2);
	int nn = nodos.size();
	for(int i=0;i<nn-1;i++){
		double dx=0, dy=0;
		
		dx = nodos[i]->x - nodos[i+1]->x;
		dy = nodos[i]->y - nodos[i+1]->y;
		
		Nodo mid = (*nodos[i] + *nodos[i+1]) / 2.0;
		
		if(dy != 0){
			M(i,0) = dx / dy;
			M(i,1) = 1.0;
			F(i) = mid.y + dx / dy * mid.x;
		}else{
			M(i,0) = 1.0;
			M(i,1) = 0;
			F(i) = mid.x;
		}
	}
	M.gauss(F,v);
	
	midPoint.x = v(0); midPoint.y = v(1);
	
	for(int i=0 ; i < nn; i++){
		Nodo mid = ( *nodos[i] + *nodos[ (i+1) % nn ] ) / 2.0;
		if( ( mid - *nodos[ (i+2) % nn ] ) * ( mid - midPoint )  < 0.001 ){
			Nodo aux = *nodos[0];
			double nsize = n.size();
			for(int i=1;i<nsize;i++)
				aux= aux + *nodos[i];
			midPoint = aux / nsize;
			break;
		}
	}
	
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
