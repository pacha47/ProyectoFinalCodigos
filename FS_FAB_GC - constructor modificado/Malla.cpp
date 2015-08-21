#include "Malla.h"
#include <cmath>
#include <stdio.h>

extern double Re, dt;

///					FUNCIONES N-S FRACTIONAL STEP
void Malla::operador_CD( double dt ){
	int n = elementos.size();
	for(int i = 0; i < n ; i++) elementos[i].operador_CD(dt);
}

void Malla::operador_CD_AB(double dt){
	int n = elementos.size();
	for(int i = 0; i < n ; i++) elementos[i].operador_CD_AB(dt);
}

void Malla::operador_P1(double dt){
	int n = elementos.size();
	double f;
	
	for(int i = 0; i < n ; i++){
		M.setRow(i, elementos[i].operador_P1(n, dt, f)); F(i) = f;}
	
	M.gradConjugado(F, P);
	
	for(int i = 0; i<n ; i++) elementos[i].setp(P(i));
}

double Malla::operador_P2(double dt){
	int n = elementos.size();
	double error=0;
	for(int i = 0; i < n ; i++) error+=elementos[i].operador_P2(dt);
	return error;
}

void Malla::primeraIte(){
	int nEle = elementos.size();
	operador_CD(dt);
	/// operador_P1(dt); a pata
	double f=0;
	
	for(int i = 0; i < nEle ; i++){
		M.setRow(elementos[i].operador_P1(nEle, dt, f)); F(i) = f;}
	
	M.gradConjugado(F, P);
	for(int i = 0; i<nEle ; i++) elementos[i].setp(P(i));
	
	operador_P2(dt);
	addvel();
	defVel();
	
}

void Malla::addvel(){ for(int i=0;i<elementos.size(); i++) elementos[i].setVelNodos() ; }
void Malla::defVel(){
	for(int i=0;i<nodos.size();i++) nodos[i].defuvp();
	
	min = max = sqrt(nodos[0].u * nodos[0].u  + nodos[0].v * nodos[0].v);
	for(int i=1;i<nodos.size();i++){
		double uv = sqrt(nodos[i].u * nodos[i].u  + nodos[i].v * nodos[i].v);
		max = (uv>max)? uv : max;
		min = (uv<min)? uv : min;
	}	
}

double Malla::iterar(double dt){
	
	operador_CD_AB(dt);
	operador_P1(dt);
	double error = operador_P2(dt);
	addvel();
	defVel();
	
	return error;
}

/// CONSTRUCTOR VACIO DE LA MALLA
Malla::Malla(){}

/// Funciones globales

bool ari_NoExiste(Nodo n1,Nodo n2){
	if(n1.aris.size() == 0 || n1.aris.size() == 0) return true;
	for(int i=0;i<n1.aris.size();i++){
		for(int j=0;j<n2.aris.size();j++){
			if(n1.aris[i] == n2.aris[j]) return false;}
	}
	return true;
}

/// CONSTRUCTOR RARO DE LA MALLA (FUNCIONA)
double Malla::makeMalla(mdouble nods, mint e, mdouble cond){
	cout<<"INICIO CARGA MALLA"<<endl;
	
	/// INICIALIZAMOS LOS NODOS 
	vector<vector<Nodo>::iterator > nodoIte;
	int nNodos = nods.size();
	for(int i=0;i<nNodos;i++)
		nodos.push_back(Nodo(nods[i][0],nods[i][1]));
	/// CON SU RESPECTIVO ITERADOR
	for(int i=0;i<nNodos;i++)
		nodoIte.push_back(nodos.begin() + i);
	cout<<"CARGO "<< nNodos<<" NODOS"<<endl;
	
	/// VARIABLE DE CONTROL (TRIANGULOS CUADRADOS)
	int nEle = e.size(), cnod = e[0].size(), id_ari = 0;
	
	/// CREAMOS LOS NUEVOS ELEMENTOS
	for(int i=0;i<nEle;i++){
		/// ASIGNAMOS LOS NODOS DEL ELEMENTO I
		vector<vector<Nodo>::iterator > nele;
		for(int n=0;n<cnod;n++)
			nele.push_back(nodoIte[e[i][n]]);
		
		/// AGREGAMOS EL NUEVO Elemento
		elementos.push_back(Elemento(nele,i));
		
		/// PARA CADA PAR DE NODOS CREAMOS UNA ARISTA Y VERIFICAMOS
		/// QUE NO HAYA SIDO CREADA Y LA GUARDAMOS EN EL ITERADOR "a"
		
		for(int j=0;j<cnod;j++){
			vector<vector<Nodo>::iterator > naux;
			/// GUARDAMOS LOS ITERADORES A LOS NODOS
			naux.push_back(nodoIte[e[i][j]]); naux.push_back(nodoIte[e[i][(j+1)%cnod]]);
			/// CREAMOS LA ARISTA CON EL I PAR DE NODOS
			Arista ari;
			ari.addNodos(naux[0],naux[1]);
			if( ari_NoExiste(*naux[0],*naux[1])){
				/// AGREGAMOS LA CONDICION DE CONTORNO DE LA Arista
				/// DEFINIMOS VARIABLES A UTILIZAR
				int ncond = cond.size();
				/// BUSCAMOS LOS NODOS EN LAS CONDICIONES "cond"
				for(int i=0;i<ncond;i++){
					if((nodos[cond[i][0]] == *naux[0] && nodos[cond[i][1]] == *naux[1]) || 
						(nodos[cond[i][0]] == *naux[1] && nodos[cond[i][1]] == *naux[0])){
							if(cond[i][2] == 1) ari.setFront(cond[i][3],cond[i][4]);
							else ari.setFront(cond[i][3]);
						}
				}
				ari.setId(id_ari);
				aristas.push_back(ari);
				naux[0]->push_ariId(id_ari);
				naux[1]->push_ariId(id_ari);
				id_ari++;
			}
		}
	}
	cout<<"CARGO ARISTAS Y ELEMENTOS"<<endl;
	nEle = elementos.size();
	int nAri = aristas.size();
	
	/// RECORREMOS LOS ELEMENTOS PARA ASIGNARLE LAS ARISTAS
	for(int i=0;i<nEle;i++){
		/// SETEAMOS LAS ARISTAS DEL ELEMENTO
		elementos[i].setAristas(aristas.begin(), elementos.begin());
	}
	cout<<"ASIGNO ARISTAS A ELEMENTOS"<<endl;
	
	for(int i=0;i<nNodos;i++)
		nodos[i].aris.clear();
	
	cout<<"ELIMINAMOS LAS ARI DE LOS NODOS"<<endl;
	
	h=aristas[0].getModulo();
	double ch = 0;
	for(int i=1;i<nAri;i++){
		h+=aristas[i].getModulo();
		ch++;
	}
	h = h/ch;
	
	
	M.setSizeCeros(nEle,cnod+2);
	F = vec(nEle); P = vec(nEle);
	
	cout<<"FIN -  CARGA MALLA"<<endl;
	
	return h;
}

void Malla::write(FILE *fs){
	fprintf(fs,"GiD Post Result File 1.0 \n");
	fprintf(fs,"Result Velocidad Fluido 1 Vector OnNodes \n");
	fprintf(fs,"Values \n");
	
	int n= nodos.size();
	for(int i=0;i<n;i++){
		fprintf(fs,"%d %.7f %.7f \n",i+1,nodos[i].u,nodos[i].v);
	}
	fprintf(fs,"End Values \n");
	
	/// RESULTADOS DE LA PRESION
	fprintf(fs,"\n");
	fprintf(fs,"Result Presion Fluido 1 Scalar OnNodes \n");
	fprintf(fs,"Values \n");
	
	for(int i=0;i<n;i++){
		fprintf(fs,"%d %.7f \n",i+1,nodos[i].p);
	}
	fprintf(fs,"End Values \n");
	
}

