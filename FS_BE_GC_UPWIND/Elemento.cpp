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
	
//	fv-= 9.8 * this->area;
	
	for(int i = 0; i < n ; i++){
		/// OBTENEMOS EL VECINO PARA SETEAR uj
		vecino = aristas[i]->getVecino(*this, dist);
		switch(aristas[i]->isFront()){
		case 1:
			uj = aristas[i]->getuv(); u_f = uj;
			/// TERMINO DIFUSIVO QUE PASA A LA FUENTE
			fu += uj.x * aristas[i]->d_ij() / Re;
			fv += uj.y * aristas[i]->d_ij() / Re;
			aP += aristas[i]->d_ij() / Re;
			/// TERMINO CONVECTIVO PASA FUENTE
			fu -= aristas[i]->getModulo() * (*aristas[i] * u_f ) * u_f.x;
			fv -= aristas[i]->getModulo() * (*aristas[i] * u_f ) * u_f.y;
			C = .0;
			break;
		case 2:
			uj = ui; u_f = uj;
			/// EL TERMINO DIFUSIVO EN PRESION IMPUESTA ES NULO
			D = 0.0;
			/// EL TERMINO CONVECTIVO PASA TODO AL ELEMENTO P
			aP += aristas[i]->getModulo() * (*aristas[i] * u_f );
			C = .0;
			/// HACEMOS QUE LA PRESION EN LA FRONTERA SEA NULA
			vecino.p = (this->p * -1.0) + aristas[i]->getp();
			break;
		default:
			/// CALCULAMOS LA VELOCIDAD EN LA FRONTERA INTERPOLANDO
			uj = Nodo(vecino.u,vecino.v);
			u_f = uj * dist[1] + ui * dist[0]; 
			/// CALCULAMOS EL VECTOR DIRECCION ENTRE LOS ELEMENTOS
			*aristas[i] % (vecino.midPoint - midPoint);
			/// CALCULAMOS EL TERMINO DIFUSIVO
			D = aristas[i]->d_ij() / Re;
			aP+=D;
			/// CALCULAMOS EL TERMINO CONVECTIVO QUICK
			vector<vector <double> > q = aristas[i]->getQuick();
			/// CALCULAMOS LA PARTE CONVECTIVA
			double vn = aristas[i]->getModulo() * (*aristas[i] * u_f );
			
			///   UPWIND
			if(vn>0){
				aP += vn;
				C = 0;
			}else{
				C = vn;
			}
			
//			/// QUICK
//			if(vn>0){
//				if(q[0][0] == this->numero){
//					aP += 				 vn  * q[0][1];
//					C   = 				 vn  * q[0][3];
//					ecuau(q[0][4]*2)   	-=vn * q[0][5];
//					ecuav(q[0][4]*2+1) 	-=vn * q[0][5];
//					ecuau(q[0][6]*2)   	-=vn * q[0][7];
//					ecuav(q[0][6]*2+1) 	-=vn * q[0][7];
//					fu += 				 vn  * q[0][8];
//					fv += 				 vn  * q[0][9]; }
//				if(q[1][0] == this->numero){
//					aP += 				 vn  * q[1][1];
//					C   = 				 vn  * q[1][3];
//					ecuau(q[1][4]*2) 	-=vn * q[1][5];
//					ecuav(q[1][4]*2+1) 	-=vn * q[1][5];
//					ecuau(q[1][6]*2) 	-=vn * q[1][7];
//					ecuav(q[1][6]*2+1) 	-=vn * q[1][7];
//					fu += 				 vn  * q[1][8];
//					fv += 				 vn  * q[1][9]; } }
//			if(vn<0){
//				if(q[0][0] == vecino.numero){
//					aP += 				 vn  * q[0][3];
//					C   = 				 vn  * q[0][1];
//					ecuau(q[0][4]*2) 	-=vn * q[0][5];
//					ecuav(q[0][4]*2+1) 	-=vn * q[0][5];
//					ecuau(q[0][6]*2) 	-=vn * q[0][7];
//					ecuav(q[0][6]*2+1) 	-=vn * q[0][7];
//					fu += 				 vn  * q[0][8];
//					fv += 				 vn  * q[0][9]; }
//				if(q[1][0] == vecino.numero){
//					aP += 				 vn  * q[1][3];
//					C   = 				 vn  * q[1][1];
//					ecuau(q[1][4]*2) 	-=vn * q[1][5];
//					ecuav(q[1][4]*2+1) 	-=vn * q[1][5];
//					ecuau(q[1][6]*2) 	-=vn * q[1][7];
//					ecuav(q[1][6]*2+1) 	-=vn * q[1][7];
//					fu += 				 vn  * q[1][8];
//					fv += 				 vn  * q[1][9]; } }
			
			break;
		}
		
		/// AGREGAMOS LA PRESION EN EL INSTANTE N
		Pu -= (vecino.p * dist[1] + this->p * dist[0] ) * aristas[i]->getModulo() * aristas[i]->getn().x;
		Pv -= (vecino.p * dist[1] + this->p * dist[0] ) * aristas[i]->getModulo() * aristas[i]->getn().y;
		
		ecuau(vecino.numero*2  ) += C - D;
		ecuav(vecino.numero*2+1) += C - D;
	}
	
	ecuau(this->numero*2  ) += aP + this->area / dt;
	ecuav(this->numero*2+1) += aP + this->area / dt;
	
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
//			ecuacion(numero)+= kdij ;
			break;
		case 2:
			vel = Nodo(u_12,v_12);
			ecuacion(numero)+= kdij;
			vecino.p = .0;
			break;
		default:
			vel = Nodo(vecino.u_12 * dist[1] + u_12 * dist[0], vecino.v_12 * dist[1] + v_12 * dist[0]);
			/// AGREGAMOS LOS VALORES CALCULADOS A LA ECUACION
			ecuacion(numero)+= kdij ;// aristas[i]->getModulo() ;
			ecuacion(vecino.numero) -= kdij ;// aristas[i]->getModulo();
			break;
		}
		*aristas[i] % (vecino.midPoint - midPoint);
		/// AGREGAMOS LA VELOCIDAD U^(1 + 1/2)
		f-= (*aristas[i] * vel) * aristas[i]->getModulo();
		
		/// AGREGAMOS LA PRESION EN EL INSTANTE ANTERIOR
		f-= (vecino.p - this->p) * kdij ;
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
			gradP *= dt / ( (*aristas[i] % dij ) * aristas[i]->getModulo());
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


vector<int> Elemento::Vecino(Elemento &P){
	vector<int> vecis;
	for(int i = 0; i<aristas.size(); i++){
		Elemento aux = aristas[i]->getVecino(*this);
		if(aux.numero != P.numero && aux.numero != this->numero)
			vecis.push_back(aux.numero);
	}
	return vecis;
}

vector<vec> Elemento::SetVecinos(int ne){
	vec vu(ne*2), vv(ne*2);
	for(int i=0;i<aristas.size();i++){
		Elemento aux = aristas[i]->getVecino(*this);
		if(aux.numero != this->numero){
			vector<int> a = aux.Vecino(*this);
			for(int j=0;j<a.size();j++)
				vu(a[j]*2) = vv(a[j]*2+1) = 1; 
			vu(aux.numero*2) = vv(aux.numero*2 + 1) = 1;
		}
	}
	vu(this->numero*2) = vv(this->numero*2 + 1) = 1;
	vector<vec> vecis; vecis.push_back(vu); vecis.push_back(vv);
	return vecis;
}

vec Elemento::SetVecinosP(int ne){
	vec vecis(ne);
	for(int i=0;i<aristas.size();i++)
		vecis((aristas[i]->getVecino(*this)).numero) = 1;
	vecis(this->numero) = 1;
	return vecis;
}


void Elemento::SetQuick(){
	Elemento N1, N2 , E;
	int na = aristas.size();
	double m1,m2,b1,b2, a1, a2, my1, my2;
	for(int i=0;i<na;i++){
		if( !aristas[i]->isFront() ){
			E = aristas[i]->getVecino(*this);
			
			/// CALCULAMOS LA PENDIENTE Y ORDENADA DE ORIGEN
			if( (E.midPoint.x - this->midPoint.x) * (E.midPoint.x - this->midPoint.x) > 1e-26 ){
					my1 = 1.0;
					m1 = (E.midPoint.y - this->midPoint.y) / (E.midPoint.x - this->midPoint.x);
					b1 = this->midPoint.y - m1 * this->midPoint.x;
			}else{	my1 = .0; m1 = -1.0;
					b1 = this->midPoint.x; }
			/// OBTENEMOS LOS ELEMENTOS VECINOS DEL W QUICK
			if( !aristas[ (i+1)%na ]->isFront() ){
					N1 = aristas[ (i+1)%na ]->getVecino(*this);
			}else{	N1 = *this;
					N1.midPoint = aristas[(i+1)%na]->getMidP() + (aristas[(i+1)%na]->getMidP() - this->midPoint); }
			if( !aristas[ (i+2)%na ]->isFront() ){ 
					N2 = aristas[ (i+2)%na ]->getVecino(*this);
			}else{	N2 = *this;
					N2.midPoint = aristas[(i+2)%na]->getMidP() + (aristas[(i+2)%na]->getMidP() - this->midPoint); }
			
			/// CALCULAMOS LA PENDIENTE Y ORDENADA DE ORIGEN
			if( ( N1.midPoint.x - N2.midPoint.x) * ( N1.midPoint.x - N2.midPoint.x) > 1e-26){
				my2 = 1.0;
				m2 = (N1.midPoint.y - N2.midPoint.y) / ( N1.midPoint.x - N2.midPoint.x);
				b2 = N1.midPoint.y - m2 * N1.midPoint.x;
			}else{	my2 = 0.0;
					m2 = -1.0;
					b2 = N1.midPoint.x;	}
			
			/// INICIALIZAMOS EL SISTEMA PARA OBTENER EL PUNTO DE INTERSECCION
			mat m(2);
			vec b(2),x(2);
			m(0,0) = my1;	m(0,1) = -m1;	b(0)   = b1;
			m(1,0) = my2;	m(1,1) = -m2;	b(1)   = b2;
			
			/// RESOLVEOS
			m.gauss(b,x);
			
			/// INICIALIZAMOS EL PUNTO DE INTERSECCION
			Nodo W(x(1),x(0));
			
			/// CALCULAMOS LOS ALPHAS DE INTERPOLACION
			m(0,0) = N1.midPoint.x;	m(0,1) = N2.midPoint.x;	b(0)   = W.x;
			m(1,0) = N1.midPoint.y;	m(1,1) = N2.midPoint.y;	b(1)   = W.y;
			/// RESOLVEMOS
			m.gauss(b,x);
			
			/// INICIALIZAMOS LOS ALPHA
			a1=x(0); a2=x(1);
			
			m(0,0) = W.x; 	m(0,1) = E.midPoint.x;	b(0) = this->midPoint.x;
			m(1,0) = W.y; 	m(1,1) = E.midPoint.y; 	b(1) = this->midPoint.y;
			
			m.gauss(b,x);
			
			double 	Qa1 = ( ( 2.0 - x(0) ) * x(1) * x(1) ) / ( 1 + x(1) - x(0) ),
				Qa2 = ( ( 1.0 - x(1) ) * (1.0 - x(0) ) * (1.0 - x(0) ) ) / ( 1.0 + x(1) - x(0) );
			
//			double Qa1 = 3.0 / 8.0, Qa2 = 1.0 / 8.0;
			
//			cout<<Qa1<<" , "<<Qa2<<" , "<<1.0-Qa1+Qa2<<endl;
			
			
//			double 	Qa1 = (aristas[i]->getMidP() - W) * (aristas[i]->getMidP() - W),
//					Qa2 = (aristas[i]->getMidP() - this->midPoint) * (aristas[i]->getMidP() - this->midPoint),
//					d   = (W -this->midPoint) * (W -this->midPoint);
//			
//			Qa1 = sqrt ( Qa1 / d);
//			Qa2 = sqrt ( Qa2 / d);
			
//			cout<<Qa1<<" , "<<Qa2<<endl;
			
			/// VERIFICAMOS FRONTERA PARA AGEGARLO A LA FUENTE
			double fu = .0 , fv = .0 , Pfront = .0;
			if(N1.numero == this->numero){
				fu += aristas[ (i+1)%na ]->getuv().x * 2.0 * a1 * Qa2;
				fv += aristas[ (i+1)%na ]->getuv().y * 2.0 * a1 * Qa2;
				Pfront += Qa2 * a1;}
			if(N2.numero == this->numero){
				fu += aristas[ (i+2)%na ]->getuv().x * 2.0 * a2 * Qa2;
				fv += aristas[ (i+2)%na ]->getuv().y * 2.0 * a2 * Qa2;
				Pfront += Qa2 * a2;}
			
//			aristas[i]->setQuick(this->numero, Pfront + Qa1, E.numero, 0, N1.numero, Qa2 * a1, N2.numero, Qa2 * a2, fu , fv);
			
			aristas[i]->setQuick(this->numero, Pfront + (1.0 - Qa1 + Qa2), E.numero, Qa1, N1.numero, Qa2 * a1, N2.numero, Qa2 * a2, fu , fv);
			aristas[i]->setQNodos(this->midPoint, E.midPoint, W, N1.midPoint , N2.midPoint);
		}
	}
	
}


void Elemento::setVelNodos(){ 
	for(int i = 0 ; i < nodos.size() ; i++) nodos[i]->adduvp(u,v,p);
}

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
	
//	for(int i=0;i < nn; i++){
//		dN.push_back( sqrt( (*nodos[i] - midPoint) * (*nodos[i] - midPoint) ) ) ;}
	
	
	
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



#include <GL/glut.h>
void Elemento::dib(){
	glColor3f(0,0,0); glLineWidth(3); glPointSize(3);
	int n=nodos.size();
	glBegin(GL_POINTS);
		for(int i=0;i<n;i++) glVertex2d(nodos[i]->x,nodos[i]->y);
	glEnd();
	
	glBegin(GL_LINES);
		for(int i=0;i<n;i++){
			glVertex2d(nodos[i]->x,nodos[i]->y);
			glVertex2d(nodos[(i+1)%n]->x,nodos[(i+1)%n]->y);
		}
	glEnd();
}
void Elemento::dibquick(vector<Elemento> &E){
	
	glColor3f(1,0,0); glLineWidth(1); glPointSize(3);
	int n=aristas.size();
//	
	vector<vector <Nodo> > q;
	
	for(int i=0;i<n;i++){
		if(!aristas[i]->isFront()){
			glColor3f( max(0,-i) , max(0, i%2) , max(0,i-1));
			
			///         0  1  2   3   4
			/// Nodos   P, E, W, W1, W2.
			q = aristas[i]->getQNodos();
			
			if(q[0][0].x == this->midPoint.x && q[0][0].y == this->midPoint.y){
				glBegin(GL_LINES);
					glVertex2d(q[0][2].x,q[0][2].y);
					glVertex2d(q[0][1].x,q[0][1].y);
					
					glVertex2d(q[0][3].x,q[0][3].y);
					glVertex2d(q[0][4].x,q[0][4].y);
				glEnd();
				
				glBegin(GL_POINTS);
					glVertex2d(q[0][0].x,q[0][0].y);
					glVertex2d(q[0][1].x,q[0][1].y);
					glVertex2d(q[0][2].x,q[0][2].y);
					glVertex2d(q[0][3].x,q[0][3].y);
					glVertex2d(q[0][4].x,q[0][4].y);
				glEnd();
			}
			
			if(q[1][0].x == this->midPoint.x && q[1][0].y == this->midPoint.y){
				glBegin(GL_LINES);
				glVertex2d(q[1][2].x,q[1][2].y);
				glVertex2d(q[1][1].x,q[1][1].y);
				
				glVertex2d(q[1][3].x,q[1][3].y);
				glVertex2d(q[1][4].x,q[1][4].y);
				glEnd();
				
				glBegin(GL_POINTS);
				glVertex2d(q[1][0].x,q[1][0].y);
				glVertex2d(q[1][1].x,q[1][1].y);
				glVertex2d(q[1][2].x,q[1][2].y);
				glVertex2d(q[1][3].x,q[1][3].y);
				glVertex2d(q[1][4].x,q[1][4].y);
				glEnd();
			}
			
		}
		
	}
	
	
	
	
	
	
	
	
	
	
//	glColor3f(1,0,0); glLineWidth(1); glPointSize(3);
//	int n=aristas.size();
////	
//	vector<vector <Nodo> > q;
//	
//	for(int i=0;i<n;i++){
//		if(!aristas[i]->isFront()){
//			glColor3f( max(0,-i) , max(0, i%2) , max(0,i-1));
//			q = aristas[i]->getQNodos();
////			cout<<q.size()<<"  "<<q[0].size()<<endl;
//			if(q[0][0] == this->numero){
////	//			q.push_back(p); q.push_back(w); q.push_back(E1); q.push_back(E2);
////	//			q.push_back(E.x); q.push_back(E.y);
////				
//				Elemento	W = E[ int(q[0][1]) ],
//							E1= E[ int(q[0][2]) ],
//							E2= E[ int(q[0][3]) ];
////				
//				
//				
//				
//				glBegin(GL_LINES);
//					glVertex2d(this->midPoint.x,this->midPoint.y);
//					glVertex2d(W.midPoint.x,W.midPoint.y);
//				glEnd();
//				
//				glBegin(GL_POINTS);
//					glVertex2d(this->midPoint.x,this->midPoint.y);
//					glVertex2d(q[0][4],q[0][5]);
//					glVertex2d(E1.midPoint.x * q[0][6] + E1.midPoint.x * q[0][7], E1.midPoint.y * q[0][6] + E1.midPoint.y * q[0][7]);
//				glEnd();
//				
//				glBegin(GL_LINES);
//					glVertex2d(E1.midPoint.x,E1.midPoint.y);
//					glVertex2d(E2.midPoint.x,E2.midPoint.y);
//				glEnd();
//				cout<<"alph1: "<<q[0][6]<<", alph2: "<<q[0][7]<<endl;
////				
//			}
//			if(q[1][0] == this->numero){
//				
//				Elemento	W = E[ q[1][1] ],
//					E1= E[ q[1][2] ],
//					E2= E[ q[1][3] ];
//				
//				glBegin(GL_LINES);
//					glVertex2d(this->midPoint.x,this->midPoint.y);
//					glVertex2d(W.midPoint.x,W.midPoint.y);
//					glVertex2d(E1.midPoint.x * q[1][6] + E1.midPoint.x * q[1][7], E1.midPoint.y * q[1][6] + E1.midPoint.y * q[1][7]);
//				glEnd();
//				
//				glBegin(GL_POINTS);
//					glVertex2d(this->midPoint.x,this->midPoint.y);
//					glVertex2d(q[1][4],q[1][5]);
//				glEnd();
//				
//				glBegin(GL_LINES);
//					glVertex2d(E1.midPoint.x,E1.midPoint.y);
//					glVertex2d(E2.midPoint.x,E2.midPoint.y);
//				glEnd();
//				
//				cout<<"alph1: "<<q[1][6]<<", alph2: "<<q[1][7]<<endl;
//				
//			}
//		
//		}
//		
//	}
}












