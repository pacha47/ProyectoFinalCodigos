#ifndef OPENGL_H
#define OPENGL_H



#include <GL/glut.h>

int w=640,h=480, el=0, boton=-1;
double escala = 400.0,escala0 = 400.0, yclick=0,xclick=0, dxt=0.01, dyt=0.01;
extern Malla m;

short modifiers=0;  // ctrl, alt, shift (de GLUT)
inline short get_modifiers() {return modifiers=(short)glutGetModifiers();}

void regen() {
	
	// matriz de proyeccion
	glMatrixMode(GL_PROJECTION);  glLoadIdentity();
	double w0=(double)w/2.0/escala,h0=(double)h/2.0/escala; // semiancho y semialto en el target
	glOrtho(-w0,w0,-h0,h0,0,1);
	glTranslated(-w0,-h0,0);
	dxt = dyt = 0;
	glMatrixMode(GL_MODELVIEW); glLoadIdentity(); // matriz del modelo->view
	glutPostRedisplay(); // avisa que hay que redibujar
}

void reshape_cb (int w, int h) {
	if (w==0||h==0) return;
	glViewport(0,0,w,h);
	regen();
	glLoadIdentity ();
}

void display_cb() {
	glClear(GL_COLOR_BUFFER_BIT);
//	double w0=(double)w/2.0/escala,h0=(double)h/2.0/escala; // semiancho y semialto en el target
//	glOrtho(-w0+dxt/2.0,w0+dxt/2.0,-h0+dyt/2.0,h0/dyt/2.0,0,1);
//	glTranslated(-w0+dxt,-h0+dyt,0);
	glTranslated(dxt,dyt,0);
	dxt = dyt = 0;
	
	
	
	m.dibele();
	m.dibquick(el);
	glutSwapBuffers();
}

void Keyboard_cb(unsigned char key,int x=0,int y=0) {
	glClear(GL_COLOR_BUFFER_BIT);
	switch(key){
	case '+':
		escala+=2;
		cout<<escala<<endl;
		break;
	case '-':
		escala-=2;
		break;
	case 'a':
		++el;
		cout<<"dib "<<el<<endl;
		
		break;
	case 'A':
		--el;
		cout<<"dib "<<el<<endl;
		break;
	case 'r':
		glLoadIdentity();
		escala=400;
		dxt = 0; dyt = 0;
		regen();
		break;
	case 'q':
		cin>>x;
		el = x;
		cout<<"dib "<<x<<endl;
		break;
	case 27: // escape => exit2
		exit(EXIT_SUCCESS);
		break;
	}
//	display_cb();
	glutPostRedisplay();
	
}

void Motion_cb(int xm, int ym){ // drag
	if (boton==GLUT_LEFT_BUTTON){
		if (modifiers==GLUT_ACTIVE_SHIFT){ // cambio de escala
			escala=escala0*exp((yclick-ym)/100.0);
			regen();
		}else{
			dxt = (xm - xclick) / escala0 / 100.0;
			dyt = (yclick - ym) / escala0 / 100.0;
			glutPostRedisplay();
		}
		
	}
//	regen();
}


void Mouse_cb(int button, int state, int x, int y){
	static bool old_rota=false;
	if (button==GLUT_LEFT_BUTTON){
		if (state==GLUT_DOWN) {
			xclick=x; yclick=y;
			boton = button;
			get_modifiers();
			glutMotionFunc(Motion_cb);
			if (modifiers==GLUT_ACTIVE_SHIFT){ // cambio de escala
				escala0=escala;
				glutPostRedisplay();
			}
		}
		else if (state==GLUT_UP){
			glutMotionFunc(0); // anula el callback para los drags
			
		}
	}
//	regen();
	glutPostRedisplay();
}


	
void initialize() {
	glutInitDisplayMode (GLUT_RGBA|GLUT_DOUBLE);
	glutInitWindowSize (w,h);
	glutInitWindowPosition (100,100);
	glutCreateWindow ("Ventana OpenGL");
	glutDisplayFunc (display_cb);
	glutReshapeFunc (reshape_cb);
	glutKeyboardFunc(Keyboard_cb);
	glutMouseFunc(Mouse_cb); // botones picados
	glClearColor(1.f,1.f,1.f,1.f);
}



#endif
