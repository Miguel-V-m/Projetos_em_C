#include <math.h>
#include <time.h>
#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include "randccp.h"

#define DEBUG //da una serie de nuneros aleatorio suscesivos
//#define SEED 316753
#define D  115 // acoplamiento
#define s2 10.0 //intensidad de ruído
#define tmax 1e6
#define test 1e2
#define samples 3 // numero de muestras independientes
#define alpha 0.5  //prescrição
#define L 50
#define N L
#define dt 1.0e-5
#define var sqrt(dt*s2)
///////////////////////////////////////////////////////
double g(double x){
  return 1+x*x;
  
}
///////////////////////////////////////////////////////
//definimos el drift com condições periodicas
double F(double x[], int i){
  return -x[i]*pow(1+x[i]*x[i],2)
    - (D/2.0)*(2*x[i]-x[(i+1)%N]-x[(i-1+N)%N])
    + (alpha-0.5)*s2*g(x[i])*2*x[i];
}
/////////////////////////////////////////////////////////
double mag(double x[]){//definição do parametro de ordem
  double soma=0.0;
  int i;

  for(i=0;i<N;i++)
    soma+=x[i];

  return soma/N;
}
/////////////////////////////////////////////////////////
double magz(double x[]){ //suma de los x^2
    double suma=0.0;
    int i;

    for(i=0;i;i++)
    suma+= x[i]*x[i];

    return suma/N;
}
////////////////////////////////////////////////////////
void rk4(double x[]){
  int i;
  double xnew[N],y[N],z[N],w[N],eta[N];

  for(i=0;i<N;i++){
    eta[i] = var*ngaussian();
    y[i] = x[i] + dt *0.5* F(x,i) + 0.5*g(x[i])*eta[i];}
  
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
///////////////////////////////////////////////////////////////
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
/////////////////////////////////////////////////////////////////////////////////////////////
///inicio del programa
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int main(void){
  double x[N];
  time_t t_ini,t_fim;
  int i,j,t,tt,h=0,m=0,s=0,tempo;
  
  t_ini = time(NULL);
  FILE *f;
  f = fopen("Cal_Mag.data","w");

  if(f==NULL){
    printf("no se ha podido habrir el archivo.\n");
    exit(1);
  }
  /////////////////////
  seed = start_randomic();
  //SEED = 316753; //se estiver definido, usa sempre a mesma semente
  //////////////////////////
  fprintf(f,"#valor del ruido: %f\n",s2);
   ///////////////////////////////////////////////q
  for (j=0;j<samples;j++){
  //AQUI SE DEFINE LA CONDICION INICIAL
  for(i=0;i<N;i++)
    //x[i]=0.1*ngaussian(); //condicion inicial aleatoria
    x[i]= 1.0;               //condicion inicial fija
    //printf("condicion inicial: %f \n ",mag(x)); //puedes descomentarlo para ver la condicion inicial 
    // (imprime la primera condicion inicial)
            
  for(t=0;t<test;t++){
        rk4(x);
     }

     for(t=0;t<tmax;t++){
     
           rk4(x);

        if(t%500==0){
        fprintf(f,"%f %f\n",t*dt,mag(x)/L);
	}

     }
     
  }
  t_fim = time(NULL);
  tempo = difftime(t_fim, t_ini);
  h = (tempo/3600);
  m = (tempo%3000)/60;
  s = (tempo%3600)%60;
  //fprintf(f,"#tempo de sumulação: %f\n",tempo); // tempo en segundos
  fprintf(f,"#tempo de sumulação: %d:%d:%d\n",h,m,s); 
    
fclose(f);
   return 0;
}
