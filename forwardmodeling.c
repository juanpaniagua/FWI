#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <mpi.h>
#define  A 10		//amplitud
#define r(i,j) 		r[Nz*(i)+(j)]
#define U1(i,j)		U1[Nz*(i)+(j)] // Un-1
#define U2(i,j)		U2[Nz*(i)+(j)] //Un
#define U3(i,j)		U3[Nz*(i)+(j)] //Un+1
#define U12(i,j)		U12[Nz*(i)+(j)] // Un-1
#define U22(i,j)		U22[Nz*(i)+(j)] //Un
#define U32(i,j)		U32[Nz*(i)+(j)] //Un+1
#define seis(i,j)	seis[Nt*i+j]
#define seis2(i,j)	seis2[Nt*i+j]
#define seisf(i,j)	seisf[Nt*i+j]
#define v(i,j)		v[nz*(i)+(j)] //velocidad
#define c(i,j)		c[Nz*(i)+(j)] //velocidad
#define c2(i,j)		c2[Nz*(i)+(j)] //velocidad
#define borde	40   //borde para fronteras absorbentes
#define alpha 0.008  //Factor de absorcion


float source(float t,float fc);
float laplac2dfd(float *indata,float *outdata, int Nz, int Nx);

int main(int argc, char *argv[])
{
int i,j,n,ns,nshots,Sx, Profsour;
int nz, nx, Nt,amp,salto,indice_archivo; //Discretizacion x, z, t
float fc,dx,dz,dt;
char caract[20],velname[20];

//======================================================================
//			Inicializacion mpi
//======================================================================
int mpi_status,mpi_rank;

mpi_status=MPI_Init(&argc,&argv);
if(mpi_status!=MPI_SUCCESS){
    fprintf(stderr,"Error inicialiando mpi, Terminando");
    MPI_Abort(MPI_COMM_WORLD,mpi_status);
}
MPI_Comm_rank(MPI_COMM_WORLD,&mpi_rank);
//Para nomcambiar nada del código de abajo, igualamos el rango de mpi al indice del disparo
ns=mpi_rank;

//======================================================================
//			Leyendo Datos del Modelo
//======================================================================

printf("Leyendo Datos del modelo \n");

FILE *sintmodel;

sintmodel = fopen("model.par","r");
    fscanf(sintmodel,"%s %s \n", &caract, &velname);
    fscanf(sintmodel,"%s %d \n", &caract, &nx);
    fscanf(sintmodel,"%s %d \n", &caract, &nz);
    fscanf(sintmodel,"%s %d \n", &caract, &Nt);
    fscanf(sintmodel,"%s %f \n", &caract, &dx);
    fscanf(sintmodel,"%s %f \n", &caract, &dz);
    fscanf(sintmodel,"%s %f \n", &caract, &dt);
    fscanf(sintmodel,"%s %f \n", &caract, &fc);
    fscanf(sintmodel,"%s %d \n", &caract, &nshots);
    fscanf(sintmodel,"%s %d \n", &caract, &Profsour);
    fscanf(sintmodel,"%s %d \n", &caract, &salto);    
fclose(sintmodel);

//======================================================================
//			Ubicación de los disparos
//======================================================================
float *Lsource=calloc(nshots,sizeof(float));

printf("Leyendo Datos de los disparos \n");

FILE *disparos;
disparos = fopen("shots.par","r");
for(i=0;i<nshots;i++){
    fscanf(disparos,"%d \n", &amp);
    Lsource[i]=amp;
}

printf("\n Datos de Disparos cargados...\n");
//======================================================================
//======================================================================
int Nz=nz+2*borde;
int Nx=nx+2*borde;

float factor;
float beta=pow(dt,2)/pow(dx,2);
float *U1 = calloc(Nx*Nz,sizeof(float));//Un-1
float *U2 = calloc(Nx*Nz,sizeof(float));//Un
float *outdata = calloc(Nx*Nz,sizeof(float));//DUn
float *U3 = calloc(Nx*Nz,sizeof(float));//Un+1
float *seis = calloc(Nt*Nx,sizeof(float));//Un+1
float  *c = calloc(Nx*Nz,sizeof(float));//velocidad
float  *c2 = calloc(Nx*Nz,sizeof(float));//velocidad
float  *r = calloc(Nx*Nz,sizeof(float));
float  *v = calloc(nx*nz,sizeof(float));//velocidad
//======================================================================
float *U12 = calloc(Nx*Nz,sizeof(float));//Un-1
float *U22 = calloc(Nx*Nz,sizeof(float));//Un
float *outdata2 = calloc(Nx*Nz,sizeof(float));//DUn
float *U32 = calloc(Nx*Nz,sizeof(float));//Un+1
float *seis2 = calloc(Nt*Nx,sizeof(float));//Un+1
float *seisf = calloc(Nt*Nx,sizeof(float));//Un+1

//======================================================================
//			Inicializa r(i,j)
//			Modelo de Velocidad
//======================================================================
printf("Leyendo Modelo de Velocidad \n");

        FILE *modelo;
        modelo = fopen (velname,"rb");
	fread(v,nx*nz*sizeof(float),1,modelo);

fclose(modelo);

//========== Inserta modelo velocidad sin bordes llenos ==========
for(i=borde;i<Nx-borde;i++){
   for(j=borde;j<Nz-borde;j++){
       c(i,j)=v(i-borde,j-borde);
       // c2(i,j)=v(i-borde,j-borde);
   }
}


//========== Borde izquierdo modelo velocidad ==========
for(i=0;i<borde;i++){
   for(j=borde;j<Nz-borde;j++){
       c(i,j)=v(0,j-borde);
      //  c2(i,j)=v(0,j-borde);
   }
}
//========== Borde derecho modelo velocidad ==========

for(i=nx+borde;i<Nx;i++){
   for(j=borde;j<Nz-borde;j++){
       c(i,j)=v(nx-1,j-borde);
     //   c2(i,j)=v(nx-1,j-borde);
   }
}


//========== Borde superior modelo velocidad ==========
for(i=0;i<Nx;i++){
   for(j=0;j<borde;j++){
       c(i,j)=c(i,borde);
      //  c2(i,j)=c2(i,borde);
   }
}

//========== Borde inferior modelo velocidad ==========
for(i=0;i<Nx;i++){
   for(j=nz+borde;j<Nz;j++){
       c(i,j)=c(i,nz+borde-1);
       // c2(i,j)=c2(i,nz+borde-1);
   }
}



for(i=0;i<Nx;i++){
   for(j=0;j<Nz;j++){
        c2(i,j)=c(i,0);
   }
}


 
printf("\n Datos de Modelo de Velocidad cargados...\n");

FILE *modelovel;
modelovel=fopen("modelovel.txt","w");
for(i=0;i<Nx;i++){
for(j=0;j<Nz;j++){

fprintf(modelovel,"%d %d %f \n", i,j,c2(i,j));
  
}
}

fclose(modelovel);


//======================================================================
//======================================================================
//for(ns=0;ns<nshots;ns++){//loop de los disparos
indice_archivo=0;
printf("Disparando la fuente %d \n",(ns+1));
Sx=Lsource[ns];
r(Sx+borde,Profsour)=1; //Ubicación de la fuente



//====================================================================
//			Programa Principal
//====================================================================
for(n=0;n<Nt;n++) { //Inicia iteración en el tiempo
printf("Calculando el campo para el disparo %d tiempo %d \n",ns+1,n);
//--------------------------------------------------------------------
//Para Obtener Un+1


laplac2dfd(U2,outdata,Nz,Nx);

//float ARP=source(n*dt,fc);

for(i=4;i<Nx-4;i++){
  for(j=4;j<Nz-4;j++)
  {
/*=================== Calcula el Un+1 ======================*/
//=========================================================================================================
U3(i,j)=2*U2(i,j)-U1(i,j)+pow(c(i,j),2)*beta*outdata[Nz*(i)+(j)]+source(n*dt,fc)*r(i,j);
  }
}
//=========================================================================================================
//					Fronteras absorbentes
//=========================================================================================================

//Para frontera x=0

for (i=0;i<borde;i++) {
    for (j=0;j<Nz;j++) {
      factor=exp(-(alpha*(borde-i))*(alpha*(borde-i)));
      U3(i,j)=U3(i,j)*factor;
      U2(i,j)=U2(i,j)*factor;
      U1(i,j)=U1(i,j)*factor;      
}
}

//Para Frontera z=0

for (i=0;i<Nx;i++){
    for (j=0;j<borde;j++) {
      factor = exp(-(alpha*(borde-j))*(alpha*(borde-j)));
      U3(i,j)= U3(i,j)*factor;
      U2(i,j)= U2(i,j)*factor;
      U1(i,j)= U1(i,j)*factor;      
}
}



//Para Frontera z=NZ
for (i=0;i<Nx;i++){
    for (j=Nz-borde;j<Nz;j++) {
      factor = exp(-(alpha *(Nz-borde-j))*(alpha*(Nz-borde-j)));
      U3(i,j)= U3(i,j)*factor;
      U2(i,j)= U2(i,j)*factor;
      U1(i,j)= U1(i,j)*factor;      
}
}
//Para Frontera x=Nx


for (j=0;j<Nz;j++) {
  for (i=Nx-borde;i<Nx;i++){
      factor = exp(-(alpha *(Nx-borde-i))*(alpha*(Nx-borde-i)));
      U3(i,j)= U3(i,j)*factor;
      U2(i,j)= U2(i,j)*factor;
      U1(i,j)= U1(i,j)*factor;      
}
}




//=========================================================================================================
//			Imprimir cada campo en cada tiempo y crear el archivo
//=========================================================================================================

  if(n%salto == 0){ 
   FILE *archvis;
char  fileoutputname[25];
sprintf(fileoutputname, "forward_%02d_%05d.bin", (ns+1), indice_archivo);
        archvis = fopen(fileoutputname,"wb");
	fwrite(U3,Nx*Nz*sizeof(float),1,archvis);
fclose(archvis); 
indice_archivo +=1;	    
}






 //========================================================================================================
//			Actualiza los campos de desplazamiento
//=========================================================================================================
for(i=0;i<Nx;i++) 
{
 for(j=0;j<Nz;j++)
  {
	U1(i,j)=U2(i,j);
	U2(i,j)=U3(i,j);
	U3(i,j)=0.0;
    }
}

for(i=0;i<Nx;i++) seis(i,n)=U2(i,borde);


} //Finaliza iteracion en el tiempo


//====================================================================
//          Nueva iteracion para restar campos
//====================================================================


for(n=0;n<Nt;n++) { //Inicia iteración en el tiempo
printf("Calculando el campo para el tiempo %d \n",n);
//--------------------------------------------------------------------
//Para Obtener Un+1


laplac2dfd(U22,outdata2,Nz,Nx);


for(i=4;i<Nx-4;i++){
  for(j=4;j<Nz-4;j++)
  {
/*=================== Calcula el Un+1 ======================*/
//=========================================================================================================
U32(i,j)=2*U22(i,j)-U12(i,j)+pow(c2(i,j),2)*beta*outdata2[Nz*(i)+(j)]+source(n*dt,fc)*r(i,j);
  }
}
//=========================================================================================================
//					Fronteras absorbentes
//=========================================================================================================

//Para frontera x=0

for (i=0;i<borde;i++) {
    for (j=0;j<Nz;j++) {
      factor=exp(-(alpha*(borde-i))*(alpha*(borde-i)));
      U32(i,j)=U32(i,j)*factor;
      U22(i,j)=U22(i,j)*factor;
      U12(i,j)=U12(i,j)*factor;      
}
}

//Para Frontera z=0

for (i=0;i<Nx;i++){
    for (j=0;j<borde;j++) {
      factor = exp(-(alpha*(borde-j))*(alpha*(borde-j)));
      U32(i,j)= U32(i,j)*factor;
      U22(i,j)= U22(i,j)*factor;
      U12(i,j)= U12(i,j)*factor;      
}
}



//Para Frontera z=NZ
for (i=0;i<Nx;i++){
    for (j=Nz-borde;j<Nz;j++) {
      factor = exp(-(alpha *(Nz-borde-j))*(alpha*(Nz-borde-j)));
      U32(i,j)= U32(i,j)*factor;
      U22(i,j)= U22(i,j)*factor;
      U12(i,j)= U12(i,j)*factor;      
}
}
//Para Frontera x=Nx


for (j=0;j<Nz;j++) {
  for (i=Nx-borde;i<Nx;i++){
      factor = exp(-(alpha *(Nx-borde-i))*(alpha*(Nx-borde-i)));
      U32(i,j)= U32(i,j)*factor;
      U22(i,j)= U22(i,j)*factor;
      U12(i,j)= U12(i,j)*factor;      
}
}

for(i=0;i<Nx;i++) 
{
 for(j=0;j<Nz;j++)
  {
	U12(i,j)=U22(i,j);
	U22(i,j)=U32(i,j);
	U32(i,j)=0.0;
    }
}

for(i=0;i<Nx;i++) seis2(i,n)=U22(i,borde);

} //Finaliza iteracion en el tiempo

//=========================================================================================================
//			Imprimir campo en cada tiempo para z=0
//=========================================================================================================
FILE *sismic;
char  fileoutputname2[25];
sprintf(fileoutputname2, "sismograma_%04d.dat", (ns+1));
sismic=fopen(fileoutputname2,"w");

for(i=0;i<Nx;i++) {
  for(n=0;n<Nt;n++){
     seisf(i,n)=seis(i,n)-seis2(i,n);
    fprintf(sismic,"%d %d %f \n", i,n,seisf(i,n));
  }
}

 
fclose(sismic);


//========================================================================

for(i=0;i<Nx;i++) {
  for(n=0;n<Nt;n++){
    seis(i,n)=0.0;
    seis2(i,n)=0.0;
    seisf(i,n)=0.0;
  }
}
for(i=0;i<Nx;i++) 
{
 for(j=0;j<Nz;j++)
  {
	U1(i,j)=0.0;
	U2(i,j)=0.0;
	U3(i,j)=0.0;
        outdata[Nz*(i)+(j)]=0.0;
        r(i,j)=0.0;
	U12(i,j)=0.0;
	U22(i,j)=0.0;
	U32(i,j)=0.0;
        outdata2[Nz*(i)+(j)]=0.0;     
    }
}

//}//finaliza iteracion fuentes


//======================================================================
//			Finalizacion mpi
//======================================================================
MPI_Finalize();

//========================================================================


}



//========================================================================
//========================================================================

//========================================================================
float source(float t,float fc)
{
  float R;
 
  float ret=0.05;
  R=A*(1-2*M_PI*M_PI*fc*fc*(t-ret)*(t-ret))*exp(-1*M_PI*M_PI*fc*fc*(t-ret)*(t-ret));
//	R=A*sin(2*M_PI*fc*(t))*exp(-2*M_PI*fc*pow((t),2));
	return R;
}

//========================================================================
//========================================================================

float laplac2dfd(float *indata,float *outdata, int Nz, int Nx)
{
  
int i,j;

for(i=4;i<Nx-4;i++){
for(j=4;j<Nz-4;j++){
outdata[Nz*(i)+(j)]=(-1.0/560)*indata[Nz*(i)+(j-4)]+(8.0/315)*indata[Nz*(i)+(j-3)]+(-1.0/5)*indata[Nz*(i)+(j-2)]+(8.0/5)*indata[Nz*(i)+(j-1)]+(-205.0/72)*indata[Nz*(i)+(j)]+(8.0/5)*indata[Nz*(i)+(j+1)]+(-1.0/5)*indata[Nz*(i)+(j+2)]+(8.0/315)*indata[Nz*(i)+(j+3)]+(-1.0/560)*indata[Nz*(i)+(j+4)];
outdata[Nz*(i)+(j)]+=(-1.0/560)*indata[Nz*(i-4)+(j)]+(8.0/315)*indata[Nz*(i-3)+(j)]+(-1.0/5)*indata[Nz*(i-2)+(j)]+(8.0/5)*indata[Nz*(i-1)+(j)]+(-205.0/72)*indata[Nz*(i)+(j)]+(8.0/5)*indata[Nz*(i+1)+(j)]+(-1.0/5)*indata[Nz*(i+2)+(j)]+(8.0/315)*indata[Nz*(i+3)+(j)]+(-1.0/560)*indata[Nz*(i+4)+(j)];

}
}

}
//========================================================================
