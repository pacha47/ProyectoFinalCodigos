
#include <iostream>
#include <cmath>

#include "Nodo.h"
#include "Elemento.h"
#include "Arista.h"

using namespace std;

extern double Re, dt;

/// ***************************** NS *****************************

/// DEFINIMOS LOS TIPOS DE FRONTERA
void Arista::setFront(double u, double v){
	oldV = v; oldU = u;
	this->u = u; this->v = v; tipoFront = 1; }
void Arista::setFront(double p){
	tipoFront = 2; this->p = p; }

/// ADIVINAMOS LOS CAMPOS DE VELOCIDAD 
void Arista::adivinar(){
	u = 0; v=0;
	oldV = v; oldU = u;}

/// SETEAMOS EL ID DE LA ARISTA
void Arista::setId(int i){id = i;}

/// REALIZAMOS LA ITERACION PARA LA ARISTA
double Arista::iterar(){
	
	/// SI ES FRONTERA SOLAMENTE RETORNAMOS EL VALOR DE LA MISMA
	if(tipoFront == 1){
		return 0;}
	
	/// EVALUAMOS SI ESTAMOS EN UNA ARISTA PARA u O v
	Nodo x(1,0);
	double alpha = fabs(x * normal);
	
	Elemento e0 = *elementos[0], e1 = *elementos[0];
	
	/// SI ES FONTERA CON PRESION IMPUESTA
	/// REFERENCIAMOS LOS ELEMENTOS 
	
	if (tipoFront == 2){
		if(alpha > 0.01 ){
			
			/// LA FRONTERA DE PRESION ESTA A LA IZQUIERDA 
			/// LADO = -1 SINO LADO = 1
			int lado = (e0.left(midPoint))? -1.0 : 1.0;
			
			this->u = e0.getAri(Nodo(.0,1)).oldV - e0.getAri(Nodo(.0,-1)).oldV + e0.getAri(Nodo(lado*-1,0)).oldU;
			return pow(u - oldU, 2);
			
		}else{
			/// LA FRONTERA DE PRESION ESTA ARRIBA
			/// LADO = -1 SINO LADO = 1
			int lado = (e0.left(midPoint))? -1.0 : 1.0;
			
			this->v = e0.getAri(Nodo(-1.0,0)).oldU - e0.getAri(Nodo(1.0,0)).oldU + e0.getAri(Nodo(0,lado*-1)).oldV;
			return pow(v - oldV, 2);
		}
	}	
	
	e0 = *elementos[0];
	e1 = *elementos[1];
	
	/// VERIFICAMOS SI ES UNA ARISTA PARA v o u
	
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
		double P = e0 - e1;
		
		/// ARISTAS LATERALES
		Arista izq = e0.getAri(Nodo(-1,0));
		Arista der = e1.getAri(Nodo(1,0));
		
		double fw = (izq.oldU + this->oldU)/2.0,
			fe = (der.oldU + this->oldU)/2.0;
		/// ARISTAS SUPERIOR E INFERIOR
		Arista top1 = e0.getAri(Nodo(0,1));
		Arista top2 = e1.getAri(Nodo(0,1));
		
		Arista bot1 = e0.getAri(Nodo(0,-1));
		Arista bot2 = e1.getAri(Nodo(0,-1));
		
		double fn = (top1.oldV + top2.oldV)/2.0,
			fs = (bot1.oldV + bot2.oldV)/2.0;
		
		/// OBTENEMOS LAS VELOCIDADES DE TOP Y BOT
		Arista top, bot;
		if(top1.tipoFront == 0) top = (top1.getVecino(e0)).getAri(Nodo(1,0));
		else top = top1;
		if(bot1.tipoFront == 0) bot = (bot1.getVecino(e0)).getAri(Nodo(1,0));
		else bot=bot1;
		
		if(top1.tipoFront == 2 || bot1.tipoFront == 2){
			this->u = (izq.oldU + der.oldU)/2.0;
			return pow(u - oldU, 2);
		}
		
		/// DEFINIMOS LAS VARIABLES DEL UPWIND
		double D = 1.0/(Re*modulo);
		double ap = D*4.0 +  max(fe,0.0) - max(-fw,0.0) + max(fn,0.0) - max(-fs,0.0) + modulo / dt,
			ae = D + max(-fe,0.0), aw = D + max(fw,0.0), 
			an = D + max(-fn,0.0), as = D + max(fs,0.0);
		
		/// ASIGNAMOS EL NUEVO VALOR DE VELOCIDAD
		u = (ae * der.oldU + aw * izq.oldU + as * bot.oldU + an * top.oldU + oldU * modulo / dt + P ) / ap;
		this->ap = ap;
		
		return pow(u - oldU, 2);
		
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
		double P = e0 - e1;
		
		/// ARISTAS SUPERIOR E INFERIOR
		Arista top = e1.getAri(Nodo(0,1));
		Arista bot = e0.getAri(Nodo(0,-1));
		
		double fn = (top.oldV + this->oldV)/2.0,
			fs = (bot.oldV + this->oldV)/2.0;
		
		/// ARISTAS LATERALES	
		Arista der1 = e0.getAri(Nodo(1,0));
		Arista der2 = e1.getAri(Nodo(1,0));
		
		Arista izq1 = e0.getAri(Nodo(-1,0));
		Arista izq2 = e1.getAri(Nodo(-1,0));
		
		double fw = (izq1.oldU + izq2.oldU)/2.0,
			fe = (der1.oldU + der2.oldU)/2.0;
		
		/// OBTENEMOS LAS VELOCIDADES DE TOP Y BOT
		Arista der, izq;
		if(der1.tipoFront == 0) der = (der1.getVecino(e0)).getAri(Nodo(0,1));
		else der = der1;
		if(izq1.tipoFront == 0) izq = (izq1.getVecino(e0)).getAri(Nodo(0,1));
		else izq=izq1;
		
		if(der1.tipoFront == 2 || izq1.tipoFront == 2){
			this->v = (bot.oldV + top.oldV)/2.0;
			return pow(v - oldV, 2);
		}
		
		/// DEFINIMOS LAS VARIABLES DEL UPWIND
		double D = 1.0/(Re*modulo);
		double ap = D*4.0 + max(fe,0.0) + max(-fw,0.0) + max(fn,0.0) + max(-fs,0.0)+ modulo  / dt,
			ae = D + max(-fe,0.0), aw = D + max(fw,0.0), 
			an = D + max(-fn,0.0), as = D + max(fs,0.0);
		
		/// ASIGNAMOS EL NUEVO VALOR DE VELOCIDAD
		v = (ae * der.oldV + aw * izq.oldV + as * bot.oldV + an * top.oldV + oldV * modulo / dt + P ) / ap;
		
		this->ap = ap;
		return pow(v - oldV, 2);
	}
}

/// CORREGIMOS LOS CAMPOS DE VELOCIDAD
double Arista::corregir(){
	
	/// SI ES FRONTERA SOLAMENTE RETORNAMOS EL VALOR DE LA MISMA
	if(tipoFront == 1){
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
			/// LA ARISTA ESTA EN LA IZQUIERDA DEL ELEMENTO e
			/// ENTONCES MODIFICAMOS LA REFERENCIA DE LOS ELEMENTOS
			Elemento aux=e0;
			e0 = e1 ; e1 = aux;}
		/// REALIZAMOS LOS PASOS PARA CORREGIR LA VELOCIDAD u		
		double u_cor = (e0 < e1) *modulo / ap;
		u = oldU + u_cor;
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


/// seteamos los vecinos para armar la matriz Ceros de U.
vec Arista::setVecinos(int n){
	
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
			/// LA ARISTA ESTA EN LA IZQUIERDA DEL ELEMENTO e
			/// ENTONCES MODIFICAMOS LA REFERENCIA DE LOS ELEMENTOS
			Elemento aux=e0;
			e0 = e1 ; e1 = aux;}
		
		/// ARISTAS LATERALES
		vecinos( (e0.getAri(Nodo(-1,0))).id ) = 1;
		vecinos( (e1.getAri(Nodo( 1,0))).id ) = 1;
		
		/// ARISTAS SUPERIOR E INFERIOR
		Arista top0 = e0.getAri(Nodo(0, 1));
		Arista bot0 = e0.getAri(Nodo(0,-1));
		
		/// OBTENEMOS LAS VELOCIDADES DE TOP Y BOT
		if(top0.tipoFront == 0) vecinos( ((top0.getVecino(e0)).getAri(Nodo(1,0)) ).id) = 1;
		if(bot0.tipoFront == 0) vecinos( ((bot0.getVecino(e0)).getAri(Nodo(1,0)) ).id) = 1;
		
		vecinos(this->id) = 1;
		
		return vecinos;
	}else{
		if( e1.top(midPoint)){ 
			/// LA ARISTA ESTA EN LA ARRIBA DEL ELEMENTO e1
			/// ENTONCES MODIFICAMOS LA REFERENCIA DE LOS ELEMENTOS
			Elemento aux=e0;
			e0 = e1 ; e1 = aux;}

		/// ARISTAS SUPERIOR E INFERIOR
		vecinos( (e1.getAri(Nodo(0, 1))).id ) = 1;
		vecinos( (e0.getAri(Nodo(0,-1))).id ) = 1;
		
		/// ARISTAS LATERALES	
		Arista der0 = e0.getAri(Nodo( 1,0));
		Arista izq0 = e0.getAri(Nodo(-1,0));
		
		/// OBTENEMOS LAS VELOCIDADES DE TOP Y BOT
		if(der0.tipoFront == 0) vecinos( ((der0.getVecino(e0)).getAri(Nodo(0,1))).id ) = 1;
		if(izq0.tipoFront == 0) vecinos( ((izq0.getVecino(e0)).getAri(Nodo(0,1))).id ) = 1;
		
		vecinos(this->id) = 1;
		
		return vecinos;
	}
	
}

/// RETORNA LA ECUACION DE VELOCIDAD DE DICHA ARISTA
vec Arista::ecuaVelocidad(int n, double &f){
	/// vector ecuacion
	vec ecua(n);
	f=0;
	/// EVALUAMOS SI ESTAMOS EN UNA ARISTA PARA u O v
	Nodo x(1,0);
	double alpha = fabs(x * normal);

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

	double 	fw = (izq.u + this->u)/2.0,
			fe = (der.u + this->u)/2.0;

	/// ARISTAS SUPERIOR E INFERIOR
	Arista top1 = e0.getAri(Nodo(0,1));
	Arista top2 = e1.getAri(Nodo(0,1));

	Arista bot1 = e0.getAri(Nodo(0,-1));
	Arista bot2 = e1.getAri(Nodo(0,-1));

	double 	fn = (top1.v + top2.v)/2.0,
			fs = (bot1.v + bot2.v)/2.0;

	///            ELEMENTO IZQUIERDO
	double ap    =  1.0 / (top1.modulo * Re) + max(-fw,0.0) * modulo; // va menos en el convectivo
	ecua(izq.id) = -1.0 / (top1.modulo * Re) - max( fw,0.0) * modulo;

	///            ELEMENTO DERECHO
	ap          +=  1.0 / (top2.modulo * Re) + max( fe,0.0) * modulo;
	ecua(der.id) = -1.0 / (top2.modulo * Re) - max(-fe,0.0) * modulo;

	///            ELEMENTO ARRIBA
	if(top1.tipoFront == 0){
		Arista top = (top1.getVecino(e0)).getAri(Nodo(1,0));

		ap          +=  1.0 / (modulo * Re) + max( fn,0.0) * top1.modulo;
		ecua(top.id) = -1.0 / (modulo * Re) - max(-fn,0.0) * top1.modulo;
	}else{
		ap += 1.0 * 2.0          / (modulo * Re) + max( top1.v,0.0) * top1.modulo;
		f  += 1.0 * 2.0 * top1.u / (modulo * Re) + max(-top1.v,0.0) * top1.modulo * top1.u;
	}

	///            ELEMENTO ABAJO
	if(bot1.tipoFront == 0){
		Arista bot = (bot1.getVecino(e0)).getAri(Nodo(1,0));

		ap          +=  1.0 / (modulo * Re) + max(-fs,0.0) * bot1.modulo; // va menos en convectivo
		ecua(bot.id) = -1.0 / (modulo * Re) - max( fs,0.0) * bot1.modulo;
	}else{
		ap += 1.0 * 2.0          / (modulo * Re) + max(-bot1.v,0.0) * bot1.modulo;
		f  += 1.0 * 2.0 * bot1.u / (modulo * Re) + max( bot1.v,0.0) * bot1.modulo * bot1.u;
	}

	///            ELEMENTO P
	ap += modulo*modulo / dt;
	f  += modulo*modulo * this->u / dt;

	ecua(this->id) = ap;

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
	
	/// GRAVEDAD
//	f-= 9.8 * modulo * modulo;
	
	/// ARISTAS SUPERIOR E INFERIOR
	Arista top = e1.getAri(Nodo(0, 1));
	Arista bot = e0.getAri(Nodo(0,-1));

	double 	fn = (top.v + this->v)/2.0,
			fs = (bot.v + this->v)/2.0;

	/// ARISTAS LATERALES	
	Arista der1 = e0.getAri(Nodo(1,0));
	Arista der2 = e1.getAri(Nodo(1,0));

	Arista izq1 = e0.getAri(Nodo(-1,0));
	Arista izq2 = e1.getAri(Nodo(-1,0));

	double 	fw = (izq1.u + izq2.u)/2.0,
			fe = (der1.u + der2.u)/2.0;


	///            ELEMENTO ARRIBA
	double ap    =  1.0 / (der2.modulo * Re) + max( fn,0.0) * modulo;
	ecua(top.id) = -1.0 / (der2.modulo * Re) - max(-fn,0.0) * modulo;

	///            ELEMENTO ABAJO
	ap          +=  1.0 / (der1.modulo * Re) + max(-fs,0.0) * modulo; // menos convectivo
	ecua(bot.id) = -1.0 / (der1.modulo * Re) - max( fs,0.0) * modulo;

	///            ELEMENTO DERECHO
	if(der1.tipoFront == 0){
		Arista der = (der1.getVecino(e0)).getAri(Nodo(0,1));

		ap          +=  1.0 / (modulo * Re) + max( fe,0.0) * der1.modulo;
		ecua(der.id) = -1.0 / (modulo * Re) - max(-fe,0.0) * der1.modulo;
	}else{
		ap += 1.0 * 2.0          / (modulo * Re) + max( der1.u,0.0) * der1.modulo;
		f  += 1.0 * 2.0 * der1.v / (modulo * Re) + max(-der1.u,0.0) * der1.modulo * der1.v;
	}

	///            ELEMENTO IZQUIERDO
	if(izq1.tipoFront == 0){
		Arista izq = (izq1.getVecino(e0)).getAri(Nodo(0,1));

		ap          +=  1.0 / (modulo * Re) + max(-fw,0.0) * izq1.modulo;  //menos convectivo
		ecua(izq.id) = -1.0 / (modulo * Re) - max( fw,0.0) * izq1.modulo;
	}else{
		ap += 1.0 * 2.0          / (modulo * Re) + max(-izq1.u,0.0) * izq1.modulo;
		f  += 1.0 * 2.0 * izq1.v / (modulo * Re) + max( izq1.u,0.0) * izq1.modulo * izq1.v;
	}


	///            ELEMENTO P
	ap += modulo*modulo / dt;
	f  += modulo*modulo * this->v / dt;

	ecua(this->id) = ap;

	this->ap = ap;

	return ecua;
	}
}//*/





void Arista::asignarUV(double uv){
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

bool Arista::operator==(Arista a){
	
	bool p1 = *nodos[0] == *a.nodos[0] && *nodos[1] == *a.nodos[1],
		p2 = *nodos[1] == *a.nodos[0] && *nodos[0] == *a.nodos[1];
	return p1 || p2;
}

double Arista::operator*(Nodo N){
	return normal*N;}

void Arista::addElemento(vector<Elemento>::iterator e){ elementos.push_back(e); }

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

Elemento Arista::getVecino(Elemento e){
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






