#include "Nodo.h"
#include "Arista.h"
#include "Elemento.h"

#include "vec.h"
#include "mat.h"

#include <cmath>
#include <iostream>
using namespace std;

extern double flechas, Re, alpha_prima;

///					FUNCIONES N-S FRACTIONAL STEP
vector<vec> Elemento::operador_CD(int ne, double dt, double &fu, double &fv){
	int n = aristas.size();
	Elemento vecino;
	Nodo ui(u,v), uj, u_f, dij;
	
	vec ecuau(2*ne);
	vec ecuav(2*ne);
	double aP=0,D=0,C=0,Pu=0,Pv=0;
	fu = fv = 0;
	
	for(int i = 0; i < n ; i++){
		/// OBTENEMOS EL VECINO PARA SETEAR uj
		vecino = aristas[i]->getVecino(*this);
		switch(aristas[i]->isFront()){
		case 1:
			uj = aristas[i]->getuv(); u_f = uj;
			/// CALCULAMOS EL VECTOR DIRECCION ENTRE LOS ELEMENTOS
			dij = aristas[i]->getMidP() - midPoint;
			vecino.midPoint = aristas[i]->getMidP();
			
			fu += uj.x * aristas[i]->getModulo() / (*aristas[i] % dij * Re);
			fv += uj.y * aristas[i]->getModulo() / (*aristas[i] % dij * Re);
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
		D = aristas[i]->getModulo() / (*aristas[i] % dij * Re);
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
		Pu -= (vecino.p + this->p) * aristas[i]->getModulo() * aristas[i]->getn().x / 2.0;
		Pv -= (vecino.p + this->p) * aristas[i]->getModulo() * aristas[i]->getn().y / 2.0;
		
		ecuau(vecino.numero*2  ) = C - D;
		ecuav(vecino.numero*2+1) = C - D;
	}
	
	ecuau(this->numero*2  ) = aP + this->area / dt;
	ecuav(this->numero*2+1) = aP + this->area / dt;
	
	fu += this->u * this->area / dt + Pu;
	fv += this->v * this->area / dt + Pv;
	
	vector<vec> ecuacion;
	ecuacion.push_back(ecuau); ecuacion.push_back(ecuav); 
	this->ap = ecuau(this->numero*2 );
	return ecuacion;
}

vector<vec> Elemento::u_prima(int ne){
	vec u(ne), v(ne);
	vector<vec> uvec;
	for(int i=0;i<aristas.size();i++){
		Elemento veci = aristas[i]->getVecino(*this);
		if(this->numero != veci.numero){
			*aristas[i] % (veci.midPoint - this->midPoint);
			u(veci.numero) += .5 * aristas[i]->getn().x;
			u(numero) += .5 * aristas[i]->getn().x;
			v(veci.numero) += .5 * aristas[i]->getn().y;
			v(numero) += .5 * aristas[i]->getn().y;
		}
	}
	u = u / ap;
	v = v / ap;
	uvec.push_back(u); uvec.push_back(v);
	return uvec;
}

vec Elemento::operador_P1(int ne, double dt, double &f){
	int n = aristas.size();
	vec ecuacion(ne);
	vector<vec> u_primaV, u_prima = this->u_prima(ne);
	f=0;
	Elemento vecino;
	Nodo vel, dij;
	for(int i = 0; i < n ; i++){
		/// OBTENEMOS EL VECINO PARA SETEAR LA VELOCIDAD EN LA FRONTERA
		vecino = aristas[i]->getVecino(*this);
		switch(aristas[i]->isFront()){
		case 1:
			vel = aristas[i]->getuv();
			break;
		case 2:
			vel = Nodo(u,v);
			break;
		default:
//			vel = Nodo( (vecino.u + u) /2.0 , (vecino.v + v)/2.0);
//			/// CALCULAMOS EL VECTOR DIRECCION ENTRE LOS ELEMENTOS
//			
//			/// CALCULAMOS LA PARTE DIFUSIVA
//			kdij =  aristas[i]->getModulo() * aristas[i]->getModulo() / (2.0 * this->ap);
			
			///   ARMA LAS ECUACIONES DE LOS P' DEL VECINO Y LOS SUMA CON EL DE P
//			std::cout<<"1 "<<std::endl;
			u_primaV = vecino.u_prima(ne);
//			std::cout<<"1 "<<std::endl;
			u_primaV[0] = u_primaV[0] + u_prima[0];
			u_primaV[1] = u_primaV[1] + u_prima[1];
			
			*aristas[i] % (vecino.midPoint - this->midPoint);
//			std::cout<<"1 "<<std::endl;
			///  MULTIPLICAMOS POR LA NORMAL
			u_primaV[0] = u_primaV[0] * aristas[i]->getn().x;
			u_primaV[1] = u_primaV[1] * aristas[i]->getn().y;
			
//			std::cout<<"1 "<<std::endl;
			///  SUMAMOS PARA OBTENER LA ECUACION DE LA ARISTA
			ecuacion = ecuacion + u_primaV[0] + u_primaV[1];
			
//			/// AGREGAMOS LOS VALORES CALCULADOS A LA ECUACION
//			ecuacion(numero)+= kdij;
//			ecuacion(vecino.numero) = -kdij;
			break;
		}
		*aristas[i] % (aristas[i]->getMidP() - this->midPoint);
		/// AGREGAMOS LA VELOCIDAD U^(1 + 1/2)
		f-= (*aristas[i] * vel) * aristas[i]->getModulo();
		
	}
	return ecuacion;
}

void Elemento::setuv(double u, double v){
	this->u = u; this->v = v;
}

void Elemento::setp_prima(double p_Prima){ 
	this->p_prima = p_Prima;
	this->p += alpha_prima*p_Prima; }

double Elemento::operador_P2(double dt){
	
	int n = aristas.size();
	double u_prima =0, v_prima =0;
	
	Elemento vecino;
	for(int i = 0; i < n ; i++){
		/// OBTENEMOS EL VECINO PARA SETEAR LA VELOCIDAD EN LA FRONTERA
		if (aristas[i]->isFront() == 0){
			*aristas[i] % (aristas[i]->getMidP() - this->midPoint);
			vecino = aristas[i]->getVecino(*this);
			u_prima += (this->p_prima + vecino.p_prima) * aristas[i]->getn().x * aristas[i]->getModulo() / (this->ap * 2.0);
			v_prima += (this->p_prima + vecino.p_prima) * aristas[i]->getn().y * aristas[i]->getModulo() / (this->ap * 2.0);
			break;
		}
		
	}
	u = u - u_prima;
	v = v - v_prima;
	
	double error = u_prima*u_prima + v_prima*v_prima;
	return error;
}


void Elemento::setVelNodos(){ 
	for(int i = 0 ; i < nodos.size() ; i++) nodos[i]->adduvp(u,v,p); }

///            CONSTRUCTOR
Elemento::Elemento(vector< vector<Nodo>::iterator > n, int id){
	///ASIGNAMOS EL ID Y LOS NODOS
	nodos = n; numero = id;
	
	/// ASIGNAMOS EL VALOR DE LAS VARIABLES
	u = v = .0;
	
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
