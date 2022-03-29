#include <math.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include "randccp.h"//puedes comentarlo
#define DEBUG
//#define SEED 316753
#define D  115 // acoplamiento
//#define s2  14.0 //intensidad de ruído
float s2;
#define tmax 1e7
#define test 1e5
#define samples 2 // numero de muestras independientes
#define alpha 0.5  //prescrição
#define L 100
#define N L
#define dt 1.0e-5
#define ds 5.0 //añadi
#define s2max 16.0
#define var sqrt(dt*s2)

////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
//definimos la difucion para el ruido multiplicativo 
double g(double x){
  return 1+x*x;
}
//////////////////////////////////////////////
/////////////////////
//definimos el drift com condições periodicas
double F(double x[], int i){
  return -x[i]*pow(1+x[i]*x[i],2)
    - D/2.0*(2*x[i]-x[(i+1)%N]-x[(i-1+N)%N])
    + (alpha-0.5)*s2*g(x[i])*2*x[i];
}
///////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////
//definição do parametro de ordem
double mag(double x[]){
  double soma=0.0;
  int i;

  for(i=0;i<N;i++)
    soma+=x[i];

  return soma/N;
}
/////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////
//definimos el mertodo de  Runge kuta de cuarta orden para la solucion de la ecuacion estocastica
void rk4(double x[]){
  int i;
  double xnew[N],y[N],z[N],w[N],eta[N];

  for(i=0;i<N;i++){
    eta[i] = var*ngaussian();
    y[i] = x[i] + dt *0.5* F(x,i) + 0.5*g(x[i])*eta[i];
  }
  for(i=0;i<N;i++){
    z[i] = x[i] + dt *0.5* F(y,i) + 0.5*g(y[i])*eta[i];}
  for(i=0;i<N;i++){
    w[i] = x[i] + dt *     F(z,i) +     g(z[i])*eta[i];}

  for(i=0;i<N;i++){
    xnew[i] = x[i] + dt*(F(x,i) + 2*F(y,i) + 2*F(z,i)+ F(w,i))/6
      + (g(x[i])+ 2*g(y[i])+ 2*g(z[i])+ g(w[i]))*eta[i]/6;}

  for(i=0;i<N;i++){
    x[i]=xnew[i];}
}
////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
// ruange kuta de segundo ordem
void rk2(double x[]){
  int i;
  double xnew[N],y[N],eta[N];

  for(i=0;i<N;i++){
    eta[i] = var*ngaussian();
    y[i] = x[i] + dt * F(x,i) + g(x[i])*eta[i];
  }
  for(i=0;i<N;i++)
    xnew[i] = x[i] + dt*(F(x,i)+F(y,i))/2 + (g(x[i])+g(y[i]))*eta[i]/2;

  for(i=0;i<N;i++)
    x[i]=xnew[i];
}
///////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
////////////////////////////////////////////////
//inicia el codigo
int main(void){
   double x[N],mm,mm2,magne,magne2,mag_media,mag2_media,varnza,desvp,standarderror,Suscepti;
  int i,j,t,tt;
  FILE *f;
  f = fopen("magXruido.data","w");

  if(f==NULL){
    printf("no se ha podido habrir el archivo.\n");
    exit(1);
  }
  /////////////////////
  seed = start_randomic();
  //SEED = 316753; //se estiver definido, usa sempre a mesma semente
  /////////////////////
  
  magne = 0;
  magne2 = 0;
  mag_media = 0;
  mag2_media = 0;
  varnza=0;
  desvp =0;
  standarderror=0;
  Suscepti=0;
  
  
   //AQUI SE DEFINE LA CONDICION INICIAL 
  for(i=0;i<N;i++)
    //x[i]=0.1*ngaussian(); //condicion inicial aleatoria
    x[i]= 1.0; //condicion inicial fija
    printf("condicion inicial: %f \n ",mag(x));
    s2=1;
  do{
  
   for (j=0;j<samples;j++){
  for(t=0;t<test;t++){
        rk4(x);
     }
     tt=0;
     mm=0;
     mm2=0;
     
     for(t=0;t<tmax;t++){
           // if(t%100==0) // puedes comentar
           // fprintf(f,"%f %f\n",t*dt,mag(x)/L);//puedes comentar
            rk4(x);
        if(t%1000==0){
	  //          mgne=mag(x);
	  //mm = mm + mag(x);
	  mm = mm + fabs(mag(x));
          mm2 = mm2 + mag(x)*mag(x);
          tt = tt + 1;
        }

     }
      magne = magne + mm/tt;
      magne2 = magne2 + mm2/tt;
  }
  
   mag_media = magne/samples;
  mag2_media = magne2/samples;
 varnza = mag2_media-mag_media*mag_media;
  //desvp =sqrt(fabs(varnza));
  standarderror = sqrt(varnza/samples);
//  Suscepti = (L*varnza)/s2;
      fprintf(f,"%f %f %f\n",s2,fabs(mag_media),standarderror);
        s2=s2+ds;
}  
 while(s2<=s2max);
fclose(f);
   return 0;
}
