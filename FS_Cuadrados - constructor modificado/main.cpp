#include <fstream>
#include <vector>
#include <iterator>
#include <iostream>
#include <cmath>
#include <string.h>

using namespace std;

typedef vector<vector<double> > mdouble;
typedef vector<vector<int> > mint;


#include <algorithm>

#include "Malla.h"

/// Variables globales
mdouble nodos, ncond;
mint ele;
Malla m;
/// VARIABLES GLOBALES :
/// iteraciones son las que hace hasta la convergencia
/// ite_vel son las que hace para estabilizar la velocidad (invento mio)
int iteraciones, ite_vel = 5;
/// max_error es el error maximo permitido entre las iteraciones
/// error_vel es un corte para la iteraciones en la velocidad (invento mio)
double Re, max_error = 1e-8, error_vel = 1e-8, dt = .0;

/// FUNCIONES ADICIONALES
void load(ifstream &f,mint &elementos, mdouble &nodos, mdouble &contorno);

int main (int argc, char **argv) {
	
	/// CARGAMOS LOS ARCHIVOS
//	ifstream file(argv[1], ios::in);
	ifstream file("DC.dat", ios::in);
	
//	FILE *fs = fopen(argv[2],"w");
	FILE *fs = fopen("DC.post.res","w");
	/// CARGAMOS LOS ELEMENTOS Y NODOS
	load(file,ele, nodos, ncond);
	
	/// ASIGNAMOS LA MALLA LEIDA
	double h = m.makeMalla(nodos,ele, ncond);
	
	dt = 0.01;
	
	cout<<endl<<"COMIENZA CALCULO CON:"<<endl;
	
	cout<<"dt : "<<dt<<", Re : "<<Re<<endl;
	cout<<"Iteraciones máximas : "<<iteraciones<<", tolerancia : "<<max_error<<endl<<endl;
	
	int ite=0;
	do{
		m.iterar();
		cout<<"Error "<<++ite<<": "<<m.error<<endl;
	}while(m.error>max_error && ite < iteraciones);
	
	m.asignar();
	
	m.write(fs);
	
//	cin>>ite;
	
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
	max_error = atof(t); /// TOLERANCIA
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
