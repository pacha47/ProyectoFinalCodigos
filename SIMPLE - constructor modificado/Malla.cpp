//
#include "Malla.h"
#include <cmath>

string s;
extern double u, alpha_p, error_vel;
extern int ite_vel, problema; 

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
//		cout<<"CARGANDO "<<i+1<<" DE "<<nEle<<" ELEMENTOS"<<endl;
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
			if(ari_NoExiste(*naux[0],*naux[1])){
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
	
	nEle = elementos.size();
	int nAri = aristas.size();
	
	cout<<"CARGO "<<nAri<<" ARISTAS Y "<<nEle<<" ELEMENTOS"<<endl;
	
	/// RECORREMOS LOS ELEMENTOS PARA ASIGNARLE LAS ARISTAS
	for(int i=0;i<nEle;i++){
		/// SETEAMOS LAS ARISTAS DEL ELEMENTO
		elementos[i].setAristas(aristas.begin(), elementos.begin());
	}
	cout<<"ASIGNO ARISTAS A ELEMENTOS"<<endl;
	
	for(int i=0;i<nNodos;i++)
		nodos[i].aris.clear();
	
	cout<<"ELIMINAMOS LAS ARI DE LOS NODOS"<<endl;
	
	h=aristas[0].modulo;;
	for(int i=1;i<nAri;i++){
		if(h>aristas[i].modulo)
			h=aristas[i].modulo;
	}
	
	cout<<"CALCULO EL H"<<endl;
	
	/// VARIABLES PARA CALCULAR EL CAMPO DE VELOCIDAD
	
	B = vec(id_ari);
	U.setSizeCeros(id_ari, cnod+2);
	
	for(int i=0;i<id_ari;i++) U.setRow(aristas[i].setVecinos(id_ari));
	
	///          VARIABLES PARA CALCULAR LA PRESION
	F = vec(nEle);
	
	K.setSizeCeros(nEle,cnod+2);
	
	for(int i=0;i<nEle;i++) K.setRow(elementos[i].setVecinos(nEle));
	
	cout<<"FIN -  CARGA MALLA"<<endl;
	
	return h;
}
/// ITERACION DEL METODO SIMPLE UTILIZANDO UPWIND
void Malla::iterar(){
	
	int an = aristas.size(), ae = elementos.size();
	double f;
	
//	int ite = 0;
//	double e = 1;
//	while ( ite++ < ite_vel  && e > error_vel){
//		e=0;
//		for(int i=0;i<an;i++)
//			e += aristas[i].iterar();
//		for(int i=0;i<an;i++)
//			aristas[i].corregir(.0);
//	}
	
	for(int i=0;i<an;i++){
		U.setRow(i,aristas[i].ecuaVelocidad(an,f));
		B(i) = f;
	}
	
	vec p_prima(an);
	U.GausSeidel(B,p_prima);
	
	for(int i=0;i<an;i++) aristas[i].asignarUV( p_prima(i) );
	
	
	for(int i=0;i<ae;i++){
		K.setRow(i, elementos[i].ecuaPresion(ae,f));
		F(i) = f;
	}
	corregir();
}

/// CORREGIMOS LAS VARIABLES DE PRESION Y VELOCIDAD
void Malla::corregir(){
	
	int an = aristas.size(), ae = elementos.size();
	vec p_prima(ae);
	
	K.gradConjugado(F,p_prima);
//	K.GausSeidel(F,p_prima);
	
	for(int i=0;i<ae;i++){
		double p_cor = p_prima(i);
		elementos[i].corregir(alpha_p, p_cor);
	}
	error = 0;
	
	for(int i=0;i<an;i++){
		error += aristas[i].corregir(1.0);}
//	asignar();
}

void Malla::asignar(){
	
	for(int i=0;i<aristas.size();i++)
		aristas[i].addvel();
	
	for(int i=0;i<elementos.size();i++)
		elementos[i].addp();
	
	for(int i=0;i<nodos.size();i++){
		nodos[i].norvel();
		nodos[i].norp();
	}
}

#include <fstream>
/// ESCRIBIMOS 
void Malla::write(FILE *fs){
	fprintf(fs,"GiD Post Result File 1.0 \n");
	fprintf(fs,"Result Velocidad Fluido 1 Vector OnNodes \n");
	fprintf(fs,"Values \n");
	
	int n= nodos.size();
	for(int i=0;i<n;i++){
		fprintf(fs,"%d %.7f %.7f \n",i+1,nodos[i].u,nodos[i].v);
	}
	fprintf(fs,"End Values \n");
	
	fprintf(fs," \n");
	
	fprintf(fs,"Result Presion Fluido 1 Scalar OnNodes \n");
	fprintf(fs,"Values \n");
	
	for(int i=0;i<n;i++){
		fprintf(fs,"%d %.7f \n",i+1,nodos[i].p);
	}
	fprintf(fs,"End Values \n");
	
}




