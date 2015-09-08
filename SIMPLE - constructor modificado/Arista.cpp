
#include <iostream>
#include <cmath>

#include "Nodo.h"
#include "Elemento.h"
#include "Arista.h"

using namespace std;

extern double Re, alpha_p , alpha_v, dt;

/// ***************************** NS *****************************

/// DEFINIMOS LOS TIPOS DE FRONTERA
void Arista::setFront(double &u, double &v){
	oldV = v; oldU = u;
	this->u = u; this->v = v; tipoFront = 1; }
void Arista::setFront(double &p){
	tipoFront = 2; this->p = p; }

/// ADIVINAMOS LOS CAMPOS DE VELOCIDAD 
void Arista::adivinar(){ u = oldU = oldV = v = .0;}

/// SETEAMOS EL ID DE LA ARISTA
void Arista::setId(int &i){id = i;}


/// seteamos los vecinos para armar la matriz Ceros de U.
vec Arista::setVecinos(int &n){
	
	vec vecinos(n);
	
	if(tipoFront == 1){
		vecinos(this->id) = 1;
		return vecinos;}
	
	/// EVALUAMOS SI ESTAMOS EN UNA ARISTA PARA u O v
	Nodo x(1,0);
	double alpha = fabs(x * normal);
	
	Elemento e0 = *elementos[0], e1 = *elementos[1];
	
	if(alpha > 0.01 ){
		if( e0.left(midPoint)){
			/// LA ARISTA ESTA EN LA IZQUIERDA DEL ELEMENTO e0
			/// ENTONCES MODIFICAMOS LA REFERENCIA DE LOS ELEMENTOS
			Elemento aux=e0;
			e0 = e1 ; e1 = aux;}
		
		/// ARISTAS LATERALES
		Arista 	W = e0.getAri(Nodo(-1,0)),
				E = e1.getAri(Nodo( 1,0));
		vecinos( W.id ) = 1;
		vecinos( E.id ) = 1;
		
		if(!W.isFront()){
			Arista WW = (W.getVecino(e0)).getAri(Nodo(-1,0));
			vecinos(WW.id) = 1;
		}
		if(!E.isFront()){
			Arista EE = (E.getVecino(e1)).getAri(Nodo(1,0));
			vecinos(EE.id) = 1;
		}
		
		/// ARISTAS SUPERIOR E INFERIOR
		Arista top0 = e0.getAri(Nodo(0, 1));
		Arista bot0 = e0.getAri(Nodo(0,-1));
		
		/// OBTENEMOS LAS VELOCIDADES DE TOP Y BOT
		if(top0.tipoFront == 0){
			Elemento TOP= top0.getVecino(e0);
			Arista 	T 	= TOP.getAri(Nodo(1,0)),
					TT	= TOP.getAri(Nodo(0,1));
			vecinos(T.id) = 1;
			if(!TT.isFront()){
				vecinos( ((TT.getVecino(TOP)).getAri(Nodo(1,0)) ).id ) = 1;
			}
		}
		if(bot0.tipoFront == 0){
			Elemento BOT= bot0.getVecino(e0);
			Arista 	B	= BOT.getAri(Nodo(1,0)),
					BB	= BOT.getAri(Nodo(0,-1));
			vecinos(B.id) = 1;
			if(!BB.isFront()){
				vecinos( ((BB.getVecino(BOT)).getAri(Nodo(1,0))).id ) = 1;
			}
		}
		
		vecinos(this->id) = 1;
		
		return vecinos;
	}else{
		if( e1.top(midPoint)){ 
			/// LA ARISTA ESTA EN LA ARRIBA DEL ELEMENTO e1
			/// ENTONCES MODIFICAMOS LA REFERENCIA DE LOS ELEMENTOS
			Elemento aux=e0;
			e0 = e1 ; e1 = aux;}

		/// ARISTAS SUPERIOR E INFERIOR
		Arista 	N = e1.getAri(Nodo(0, 1)),
				S = e0.getAri(Nodo(0,-1));
		vecinos( N.id ) = 1;
		vecinos( S.id ) = 1;
		
		if(!N.isFront()){
			Arista NN = (N.getVecino(e1)).getAri(Nodo(0,1));
			vecinos(NN.id) = 1;
		}
		if(!S.isFront()){
			Arista SS = (S.getVecino(e0)).getAri(Nodo(0,-1));
			vecinos(SS.id) = 1;
		}
		
		/// ARISTAS LATERALES	
		Arista der0 = e0.getAri(Nodo( 1,0));
		Arista izq0 = e0.getAri(Nodo(-1,0));
		
		/// OBTENEMOS LAS VELOCIDADES DE TOP Y BOT
		if(der0.tipoFront == 0){
			Elemento DER= der0.getVecino(e0);
			Arista D = DER.getAri(Nodo(0,1)),
				DD = DER.getAri(Nodo(1,0));
			vecinos (D.id) = 1;
			if(!DD.isFront()){
				vecinos( ((DD.getVecino(DER)).getAri(Nodo(0,1))).id) = 1;
			}
		}
		if(izq0.tipoFront == 0){
			Elemento IZQ= izq0.getVecino(e0);
			Arista 	I = IZQ.getAri(Nodo( 0,1)),
					II= IZQ.getAri(Nodo(-1,0));
			vecinos(I.id) = 1;
			if(!II.isFront()){
				vecinos( ((II.getVecino(IZQ)).getAri(Nodo(0,1)) ).id ) = 1;
			}
			vecinos( ((izq0.getVecino(e0)).getAri(Nodo(0,1))).id ) = 1;
		}
		
		vecinos(this->id) = 1;
		
		return vecinos;
	}
	
}

/// RETORNA LA ECUACION DE VELOCIDAD DE DICHA ARISTA
vec Arista::ecuaVelocidad(int &n, double &f){
	/// vector ecuacion
	vec ecua(n);
	f=0;
	/// EVALUAMOS SI ESTAMOS EN UNA ARISTA PARA u O v
	Nodo x(1,0);
	double alpha = fabs(x * normal), 
		a1 = 3.0/8.0, a2 = 1.0/8.0, a3 = 6.0/8.0;

	if(tipoFront == 1){
		ecua(this->id) = 1;
		f = (alpha > 0.01) ? u : v;
		this->ap = 1;
		return ecua;}

	Elemento e0 = *elementos[0], e1 = *elementos[1];

	if(alpha > 0.01 ){
	if( e0.left(midPoint)){
		/// LA ARISTA ESTA EN LA IZQUIERDA DEL ELEMENTO e
		/// ENTONCES MODIFICAMOS LA REFERENCIA DE LOS ELEMENTOS
		Elemento aux=e0;
		e0 = e1 ; e1 = aux;}
	/// REALIZAMOS LOS PASOS PARA ADQUIRIR LOS COEFICIENTES
	/// NECESARIOS PARA REALIZAR LA ITERACION

	///  PARA LAS VELOCIDADES EN u EL ELEMENTO e0 ESTA A LA IZQUIERDA
	/// Y EL ELEMENTO e1 ESTA A LA DERECHA SIN IMPORTAR SI ES FRONTERA O NO

	/// GRADIENTE DE PRESION
	f = (e0 - e1) * modulo;
	
	/// ARISTAS LATERALES
	Arista izq = e0.getAri(Nodo(-1,0));
	Arista der = e1.getAri(Nodo( 1,0));

	

	/// ARISTAS SUPERIOR E INFERIOR
	Arista top1 = e0.getAri(Nodo(0,1));
	Arista top2 = e1.getAri(Nodo(0,1));

	Arista bot1 = e0.getAri(Nodo(0,-1));
	Arista bot2 = e1.getAri(Nodo(0,-1));

	
	///			FLUJOS DIFUSIVOS
	
	///            ELEMENTO IZQUIERDO
	double ap    =  1.0 / Re;
	ecua(izq.id) = -1.0 / Re;
	
	///            ELEMENTO DERECHO
	ap          +=  1.0 / Re;
	ecua(der.id) = -1.0 / Re;

	///            ELEMENTO ARRIBA
	if(top1.tipoFront == 0){
		Arista top = (top1.getVecino(e0)).getAri(Nodo(1,0));
		ap          +=  1.0 / Re;
		ecua(top.id) = -1.0 / Re;
	}else{
		ap += 2.0 	/ Re;
		f  += 2.0 * top1.u / Re;
	}

	///            ELEMENTO ABAJO
	if(bot1.tipoFront == 0){
		Arista bot = (bot1.getVecino(e0)).getAri(Nodo(1,0));
		ap          +=  1.0 / Re;
		ecua(bot.id) = -1.0 / Re;
	}else{
		ap += 2.0 	/ Re;
		f  += 2.0 * bot1.u / Re;
	}
	
	/// 				FLUJOS CONVECTIVOS
	double 	fw = (izq.u + this->u)/2.0 * -1.0 * modulo,
			fe = (der.u + this->u)/2.0 *  1.0 * modulo,
			fn = (top1.v + top2.v)/2.0 *  1.0 * modulo,
			fs = (bot1.v + bot2.v)/2.0 * -1.0 * modulo;
	
	/// 					CARA W
	if ( fw > 0 ){
		ecua(der.id)-= a2 * fw ;
		ecua(izq.id)+= a1 * fw ;
		ap			+= a3 * fw ;
	}else{
		ecua(izq.id)+= a3 * fw;
		ap			+= a1 * fw;
		if(!izq.isFront()){
			Arista WW = izq.getVecino(e0).getAri(Nodo(-1,0));
			ecua(WW.id) = -a2 * fw;
		}else{
			ecua(izq.id)+= a2 * fw;
			f			+= izq.u * 2.0 * a2 * fw;
		}
	}
	
	/// 					CARA E
	if ( fe > 0 ){
		ecua(der.id)+= a1 * fe;
		ecua(izq.id)-= a2 * fe;
		ap			+= a3 * fe;
	}else{
		ecua(der.id)+= a3 * fe;
		ap			+= a1 * fe;
		if(!der.isFront()){
			Arista EE = der.getVecino(e1).getAri(Nodo(1,0));
			ecua(EE.id) = -a2 * fe;
		}else{
			ecua(der.id) += a2 * fe;
			f += der.u * 2.0 * a2 * fe;
		}
	}
	
	Arista 	top = (top1.getVecino(e0)).getAri(Nodo(1,0)),
			bot = (bot1.getVecino(e0)).getAri(Nodo(1,0));
	/// 					CARA N
	if(!top1.isFront()){
		if ( fn > 0 ){
			ecua(top.id)+= a1 * fn;
			ap			+= a3 * fn;
			if(!bot1.isFront()){
				ecua(bot.id)-= a2 * fn;
			}else{
				f+= bot1.u * 2.0 * a2 * fn;
				ap += a2 * fn;
			}
		}else{
			ecua(top.id)+= a3 * fn;
			ap			+= a1 * fn;
			Elemento T = top1.getVecino(e0);
			Arista NN = T.getAri(Nodo(0,1));
			if(!NN.isFront()){
					NN = (NN.getVecino(T)).getAri(Nodo(1,0));
					ecua(NN.id) = -a2 * fn;
				}else{
					f+= NN.u * 2.0 * a2 * fn;
					ecua(top.id) += a2 * fn;
				}
		}
	}
	
	/// 					CARA S
	if(!bot1.isFront()){
		if ( fs > 0 ){
			ecua(bot.id)+= a1 * fs;
			ap			+= a3 * fs;
			if(!top1.isFront()){
				ecua(top.id)-= a2 * fs;
			}else{
				f+= top1.u * 2.0 * a2 * fn;
				ap+= a2 * fn;
			}
		}else{
			ecua(bot.id)+= a3 * fs;
			ap			+= a1 * fs;
			Elemento T = bot1.getVecino(e0);
			Arista SS = T.getAri(Nodo(0,-1));
			if(!SS.isFront()){
				SS = (SS.getVecino(T)).getAri(Nodo(1,0));
				ecua(SS.id) = -a2 * fs;
			}else{
				f+= SS.u * 2.0 * a2 * fn;
				ecua(bot.id) += a2 * fs;
			}
		}
	}
	
	///            ELEMENTO P
	ap += modulo*modulo / dt;
	f  += modulo*modulo * this->u / dt;
	
	/// RELAJACION
//	f+= ( 1.0 - alpha_v ) * ap * this->u / alpha_v;
	
	ecua(this->id) = ap ;// alpha_v;

	this->ap = ap;

	return ecua;

	}else{
	if( e1.top(midPoint)){ 
		/// LA ARISTA ESTA EN LA ARRIBA DEL ELEMENTO e1
		/// ENTONCES MODIFICAMOS LA REFERENCIA DE LOS ELEMENTOS
		Elemento aux=e0;
		e0 = e1 ; e1 = aux;}
	/// REALIZAMOS LOS PASOS PARA ADQUIRIR LOS COEFICIENTES
	/// NECESARIOS PARA REALIZAR LA ITERACION

	/// PARA LAS VELOCIDADES EN v EL ELEMENTO e0 ESTA ABAJO Y EL
	/// ELEMENTO e1 ESTA ARRIBA SIN IMPORTAR SI ES FRONTERA O NO

	/// GRADIENTE DE PRESION
	f = (e0 - e1) * modulo;
	
	/// ARISTAS SUPERIOR E INFERIOR
	Arista top = e1.getAri(Nodo(0, 1));
	Arista bot = e0.getAri(Nodo(0,-1));

	/// ARISTAS LATERALES	
	Arista der1 = e0.getAri(Nodo(1,0));
	Arista der2 = e1.getAri(Nodo(1,0));

	Arista izq1 = e0.getAri(Nodo(-1,0));
	Arista izq2 = e1.getAri(Nodo(-1,0));

	///            ELEMENTO ARRIBA
	double ap    =  1.0 / Re;
	ecua(top.id) = -1.0 / Re;

	///            ELEMENTO ABAJO
	ap          +=  1.0 / Re;
	ecua(bot.id) = -1.0 / Re;

	///            ELEMENTO DERECHO
	if(der1.tipoFront == 0){
		Arista der = (der1.getVecino(e0)).getAri(Nodo(0,1));

		ap          +=  1.0 / Re;
		ecua(der.id) = -1.0 / Re;
	}else{
		ap += 2.0 	/ Re;
		f  += 2.0 * der1.v / Re;
	}

	///            ELEMENTO IZQUIERDO
	if(izq1.tipoFront == 0){
		Arista izq = (izq1.getVecino(e0)).getAri(Nodo(0,1));

		ap          +=  1.0 / Re;
		ecua(izq.id) = -1.0 / Re;
	}else{
		ap += 2.0 	/ Re;
		f  += 2.0 * izq1.v / Re;
	}
	
	/// 				FLUJOS CONVECTIVOS
	double 	fn = (top.v + this->v)/2.0 *  1.0 * modulo,
			fs = (bot.v + this->v)/2.0 * -1.0 * modulo,
			fw = (izq1.u + izq2.u)/2.0 * -1.0 * modulo,
			fe = (der1.u + der2.u)/2.0 *  1.0 * modulo;
	
	/// 					CARA N
	if ( fn > 0 ){
		ecua(bot.id)-= a2 * fn;
		ecua(top.id)+= a1 * fn;
		ap			+= a3 * fn;
	}else{
		ecua(top.id)+= a3 * fn;
		ap			+= a1 * fn;
		if(!top.isFront()){
			Arista NN = top.getVecino(e1).getAri(Nodo(0,1));
			ecua(NN.id) = -a2 * fn;
		}else{
			ecua(top.id) += a2 * fn;
			f +=  top.v * 2.0 * a2 * fn;
		}
	}
	
	/// 					CARA S
	if ( fs > 0 ){
		ecua(bot.id)+= a1 * fs;
		ecua(top.id)-= a2 * fs;
		ap			+= a3 * fs;
	}else{
		ecua(bot.id)+= a3 * fs;
		ap			+= a1 * fs;
		if(!bot.isFront()){
			Arista SS = bot.getVecino(e0).getAri(Nodo(0,-1));
			ecua(SS.id) = -a2 * fs;
		}else{
			ecua(bot.id) += a2 * fs;
			f += bot.v * 2.0 * a2 * fs;
		}
	}
	
	Arista 	der = (der1.getVecino(e0)).getAri(Nodo(0,1)),
		izq = (izq1.getVecino(e0)).getAri(Nodo(0,1));
	/// 					CARA E
	if(!der1.isFront()){
		if ( fe > 0 ){
			ecua(der.id)+= a1 * fe;
			ap			+= a3 * fe;
			if(!izq1.isFront()){
				ecua(izq.id)-= a2 * fe;
			}else{
				f+= izq1.v * 2.0 * a2 * fe;
				ap+= a2 * fe;
			}
		}else{
			ecua(der.id)+= a3 * fe;
			ap			+= a1 * fe;
			Elemento T = der1.getVecino(e0);
			Arista EE = T.getAri(Nodo(1,0));
			if(!EE.isFront()){
				EE = (EE.getVecino(T)).getAri(Nodo(0,1));
				ecua(EE.id) = -a2 * fe;
			}else{
				f+= EE.v * 2.0 * a2 * fe;
				ecua(der.id) += a2 * fe;
			}
		}
	}
	
	/// 					CARA W
	if(!izq1.isFront()){
		if ( fw > 0 ){
			ecua(izq.id)+= a1 * fw;
			ap			+= a3 * fw;
			if(!der1.isFront()){
				ecua(der.id)-= a2 * fw;
			}else{
				f+= der1.v * 2.0 * a2 * fe;
				ap += a2 * fw;
			}
		}else{
			ecua(izq.id)+= a3 * fw;
			ap			+= a1 * fw;
			Elemento T = izq1.getVecino(e0);
			Arista WW = T.getAri(Nodo(-1,0));
			if(!WW.isFront()){
				WW = (WW.getVecino(T)).getAri(Nodo(0,1));
				ecua(WW.id) = -a2 * fw;
			}else{
				f+= WW.v *2.0 * a2 * fw;
				ecua(izq.id) += a2 * fw;
			}
		}
	}

	///            ELEMENTO P
	ap += modulo*modulo / dt;
	f  += modulo*modulo * this->v / dt;
	
	/// RELAJACION
//	f+= ( 1.0 - alpha_v ) * ap * this->v / alpha_v;
	
	ecua(this->id) = ap ;// alpha_v;

	this->ap = ap;

	return ecua;
	}
}//*/


/// CORREGIMOS LOS CAMPOS DE VELOCIDAD
double Arista::corregir(){
	
	/// SI ES FRONTERA SOLAMENTE RETORNAMOS EL VALOR DE LA MISMA
	if(tipoFront == 1){
//		cout<<endl;
		return .0;}
	oldV=v; oldU=u;
	
	/// EVALUAMOS SI ESTAMOS EN UNA ARISTA PARA u O v
	Nodo x(1,0);
	double xy = fabs(x * normal);
	
	/// REFERENCIAMOS LOS ELEMENTOS 
	Elemento e0 = *elementos[0], e1 = *elementos[0];
	
	/// REFERENCIAMOS LOS ELEMENTOS SI NO SON FRONTERA
	if(tipoFront == 0 ){
		e0 = *elementos[0];
		e1 = *elementos[1];
	}
	
	/// VERIFICAMOS SI ES UNA ARISTA PARA v o u
	if( xy > 0.01){
		if( e0.left(midPoint) ){ 
			/// LA ARISTA ESTA EN LA IZQUIERDA DEL ELEMENTO e0
			/// ENTONCES MODIFICAMOS LA REFERENCIA DE LOS ELEMENTOS
			Elemento aux=e0;
			e0 = e1 ; e1 = aux;}
		/// REALIZAMOS LOS PASOS PARA CORREGIR LA VELOCIDAD u		
		double u_cor = (e0 < e1) * modulo / ap;
		u = oldU + u_cor;
//		cout<<u<<endl;
		return pow(oldU - u,2.0);
	}else{
		if( e1.top(midPoint) ){ 
			/// LA ARISTA ESTA EN LA ARRIBA DEL ELEMENTO e
			/// ENTONCES MODIFICAMOS LA REFERENCIA DE LOS ELEMENTOS
			Elemento aux=e0;
			e0 = e1 ; e1 = aux;}
		/// REALIZAMOS LOS PASOS PARA CORREGIR LA VELOCIDAD v
		double v_cor = (e0 < e1) * modulo / ap ;
		v= oldV + v_cor;
//		cout<<v<<endl;
		return pow(oldV - v,2.0);
	}
}

void Arista::addvel(){
	/// EVALUAMOS SI ESTAMOS EN UNA ARISTA PARA u O v
	Nodo x(1,0);
	double alpha = x * normal;
	if(tipoFront == 1){
		double n = 100;
		nodos[0]->cv += n;
		nodos[1]->cv += n;
		nodos[0]->v += v*n;
		nodos[1]->v += v*n;
		nodos[0]->cu += n;
		nodos[1]->cu += n;
		nodos[0]->u += u*n;
		nodos[1]->u += u*n;
		return; }
	if(fabs(alpha) > 0.0001){
		nodos[0]->cu ++;
		nodos[1]->cu ++;
		nodos[0]->u += u;
		nodos[1]->u += u;
		return;}
	nodos[0]->cv ++;
	nodos[1]->cv ++;
	nodos[0]->v += v;
	nodos[1]->v += v;
}




void Arista::asignarUV(double &uv){
	if(tipoFront==1) return;
	/// EVALUAMOS SI ESTAMOS EN UNA ARISTA PARA u O v
	Nodo x(1,0);
	double alpha = fabs(x * normal);
	if(alpha > 0.01){
		this->u = uv;
	}else{
		this->v = uv;}
}



/// ************************************************************** 

bool Arista::operator==(Arista &a){
	
	bool p1 = *nodos[0] == *a.nodos[0] && *nodos[1] == *a.nodos[1],
		p2 = *nodos[1] == *a.nodos[0] && *nodos[0] == *a.nodos[1];
	return p1 || p2;
}

double Arista::operator*(Nodo &N){
	return normal*N;}

void Arista::addElemento(vector<Elemento>::iterator e){ elementos.push_back(e); }

void Arista::addNodos(vector<Nodo>::iterator &n1,vector<Nodo>::iterator &n2){
	/// ASIGNAMOS LOS NODOS
	nodos.push_back(n1); nodos.push_back(n2);
	/// CALCULAMOS LA NORMAL
	normal.x = (*n1).y - (*n2).y;
	normal.y = (*n2).x - (*n1).x;
	/// NORMALIZAMOS LA NORMAL Y ASIGNAMOS EL MODULO DE LA ARISTA 
	normal.modulo = modulo = sqrt(normal.x * normal.x + normal.y * normal.y);
//	modulo = sqrt ( (*n1 - *n2) * (*n1 - *n2) );
	normal.x = normal.x / normal.modulo;
	normal.y = normal.y / normal.modulo;
	/// RECALCULAMOS EL MODULO == 1
	normal.modulo = normal.x * normal.x + normal.y * normal.y;
	/// CALCULAMOS EL PUNTO MEDIO
	midPoint = (*nodos[0] + *nodos[1]) / 2.0;
}

Elemento Arista::getVecino(Elemento &e){
	if(tipoFront == 1)
		return e;
	if(tipoFront == 2){
		Elemento aux = e;
		e.setPresion(p);
		return aux;
	}
	if(e == (*elementos[0]))
		return (*elementos[1]);
	return (*elementos[0]);
}

double Arista::getmodulo(){
	if(modulo!=-1)
		return modulo;
	modulo = sqrt((*nodos[0] - *nodos[1])*(*nodos[0] - *nodos[1]));
	return modulo;}






