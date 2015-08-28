#include <fstream>
#include <vector>
#include <iterator>
#include <cmath>
#include <cstring>

using namespace std;

typedef vector<vector<double> > mdouble;
typedef vector<vector<int> > mint;


#include <algorithm>
#include "Malla.h"

/// Variables globales
mdouble nodos, ncond;
mint ele;
Malla m;
int iteraciones;
double Re, dt, e;

void load(ifstream &f,mint &elementos, mdouble &nodos, mdouble &contorno);

//FILE *fs;

int main (int argc, char **argv) {
	
	/// CARGAMOS LOS ARCHIVOS
//	ifstream file(argv[1], ios::in);
	ifstream file("Cuadrado crusado TRIAN.dat", ios::in);
	
//	FILE *fs = fopen(argv[2],"w");
	FILE *fs = fopen("Cuadrado crusado TRIAN.post.res","w");
	
	/// CARGAMOS LOS ELEMENTOS Y NODOS
	load(file,ele, nodos, ncond);
	
	/// ASIGNAMOS LA MALLA LEIDA
	m.makeMalla(nodos,ele, ncond);
	
	cout<<endl<<"COMIENZA CALCULO CON:"<<endl;
	
	cout<<"dt : "<<dt<<", Re : "<<Re<<endl;
	cout<<"Iteraciones máximas : "<<iteraciones<<", tolerancia : "<<e<<endl<<endl;
	
	m.primeraIte();
	
	iteraciones = 6000;
	
	int i=0, ii=2;
	double ite_e=10;
//	m.write(fs);
	while (ite_e>e && iteraciones>i){
		ite_e = m.iterar(dt);
		cout<<"Error "<<++i<<": "<<ite_e<<endl;
//		if( i%5 == 0 ) {
//			m.addvel();
//			m.defVel();
//			m.write(fs,ii++);
//		}
	}
	
	m.addvel();
	m.defVel();
	m.write(fs);
	
	cin>>i;
	
	return 0;
}


void load(ifstream &f,mint &elementos, mdouble &nodos, mdouble &contorno){
	
	nodos.clear();
	elementos.clear();
	contorno.clear();
	/// VARIABLES
	int n,o,p;
	char t[100];
	double x;
	
	/// CARGAMOS LAS VARIABLES GENERALES DEL PROBLEMA
	f.getline(t,100);
	e = atof(t); /// TOLERANCIA
	f.getline(t,100);
	iteraciones = atoi(t); /// ITERACIONES MAXIMAS
	f.getline(t,100);
	dt = atof(t); /// PASO DE TIEMPO
	f.getline(t,100);
	Re = atof(t); /// NUMERO DE REYNOLDS
	
	/// LEEMOS LA CANTIDAD DE NODOS
	f.getline(t,100);
	n = atoi(t);
	/// ASIGNAMOS LOS NODOS
	for(int i=0;i<n;i++){
		vector<double> nod;
		/// COORDENADA X
		f.getline(t,100);
		x = atof(t);
		nod.push_back(x);
		/// COORDENADA Y
		f.getline(t,100);
		x = atof(t);
		nod.push_back(x);
		nodos.push_back(nod);
	}
	
	/// LEEMOS LA CANTIDAD DE ELEMENTOS
	f.getline(t,100);
	n = atoi(t);
	
	/// LEEMOS LA CANTIDAD DE NODOS POR ELEMENTO
	f.getline(t,100);
	o = atoi(t);
	
	/// LEEMOS LOS INDICES DE LOS ELEMENTOS
	for(int i=0;i<n;i++){
		vector<int> e;
		f.getline(t,100);
		char *ptr;
		/// LOS NUMEROS DE LOS NODOS ESTAN SEPARADOS POR UN ESPACIO
		ptr = strtok(t," ");
		p = atoi(ptr);
		e.push_back(p-1);
		for(int j=0;j<o-1;j++){
			// mandarle un null parece que hace que continue con el puntero del anterior strtok
			ptr = strtok(NULL," ");
			p = atoi(ptr);
			e.push_back(p-1);
		}
		elementos.push_back(e);
	}
	/// LEEMOS LAS CONDICIONES DE CONTORNO
	
	/// LA ORGANIZACION DEL CONTORNO ES
	/// NumNodo TipoFront Valores
	
	/// LEEMOS LA CANTIDAD DE NODOS CON WALL
	f.getline(t,100);
	n = atoi(t);
	for(int i=0;i<n;i++){
		vector<double> nContorno;
		/// NODO 1
		f.getline(t,100);
		nContorno.push_back(atoi(t)-1);
		/// NODO 2
		f.getline(t,100);
		nContorno.push_back(atoi(t)-1);
		/// TIPO DE CONTORNO
		nContorno.push_back(1.0);
		/// ASIGNAMOS u = v = 0 porque es wall
		nContorno.push_back(.0);
		nContorno.push_back(.0);
		contorno.push_back(nContorno);
	}
	
	/// LEEMOS LA CANTIDAD DE NODOS CON VELOCIDAD IMPUESTA
	f.getline(t,100);
	n = atoi(t);
	for(int i=0;i<n;i++){
		vector<double> nContorno;
		/// NODO 1
		f.getline(t,100);
		nContorno.push_back(atoi(t)-1);
		/// NODO 2
		f.getline(t,100);
		nContorno.push_back(atoi(t)-1);
		/// TIPO DE CONTORNO
		nContorno.push_back(1.0);
		/// LEEMOS u
		f.getline(t,100);
		nContorno.push_back(atof(t));
		/// LEEMOS v
		f.getline(t,100);
		nContorno.push_back(atof(t));
		contorno.push_back(nContorno);
	}
	
	/// LEEMOS LA CANTIDAD DE NODOS CON PRESION IMPUESTA
	f.getline(t,100);
	n = atoi(t);
	for(int i=0;i<n;i++){
		vector<double> nContorno;
		/// NODO 1
		f.getline(t,100);
		nContorno.push_back(atof(t)-1);
		/// NODO 2
		f.getline(t,100);
		nContorno.push_back(atof(t)-1);
		/// TIPO DE CONTORNO
		nContorno.push_back(2.0);
		/// LEEMOS LA PRESION
		f.getline(t,100);
		nContorno.push_back(atof(t));
		contorno.push_back(nContorno);
	}
	
	f.close();
}
