# include<math.h>
# include<string.h>
# include<stdlib.h>
# include<conio.h>
# include<graphics.h>
# include"..\common.h"
# include"..\rparam\rparam.h"
# include"..\iwork\lab\plofunc.h"
# include"four.h"

float CheckStep=0.1;
float Amplitude,Period,Phase,Noise;
float *Data;
float XCoeff,Xmin,YCoeff,Dx,Dstep;
FILE *InFile,*OutFile=NULL;
int InType=0;

float Signal(float x){
 float y;
 static float xgrp=-1,xprev=-1,xc,yprev,yc;

 if(InType==0){
  y=Amplitude*sin(2*M_PI*x/Period+M_PI*Phase)+Noise*(16000.-random(32000))/16000.;
 }
 else{
  if(xprev<0){
   if(InType==2)fscanf(InFile,"%f %f",&xc,&yc);
   else fscanf(InFile,"%f",&yc);
   xc=xprev=0;
   yc=yprev=yc*YCoeff;
  }
  while(xc<=x){
   xprev=xc;
   yprev=yc;
   if(InType==2){
    fscanf(InFile,"%f %f",&xc,&yc);
    xc=XCoeff*(xc-Xmin);
   }
   else{
    fscanf(InFile,"%f",&yc);
    xc+=Dx;
   }
   yc=yc*YCoeff;
  }
  y=lin_approx(xprev,yprev,xc,yc,x);
 }
 Data[2*(int)(x/CheckStep)]=y;

 if(OutFile && x-xgrp>Dstep){
  fprintf(OutFile,"%f %f\n",x,y);
  xgrp=x;
 }

 return y;
}

float ResFreq=1.,ResAsin,ResAcos;
float ResSignal(float x){
 x-=CheckStep;
 return ResAcos*cos(2*M_PI*x*ResFreq)+ResAsin*sin(2*M_PI*x*ResFreq);
}

int N;
float All(float x){
 x-=CheckStep;
 float y=0;
 for(int i=0;i<=N/2;i++){
  y+=Data[2*i]*cos(2*M_PI*x*i)+Data[2*i+1]*sin(2*M_PI*x*i);
 }
 for(i=-1;i>=-N/2;i--){
  y+=Data[2*N+2*i]*cos(2*M_PI*x*i)+Data[2*N+(2*i+1)]*sin(2*M_PI*x*i);
 }
 return y/N;
}

void main(){
 long n;
 Open_param_file("filter.cfg");
 Read_param("p: %ld",&n);
 Set_stop(0);
 float s=0;
 Read_param("Сигнал/шум: %f",&s);
 Noise=1./(1.+s);
 Amplitude=Noise*s;
 Read_param("Период сигнала: %f",&Period);
 Read_param("Фазовый сдвиг: %f",&Phase);

 char in_file[50]="",out_file[50]="";

 Set_comment_char('/');
 if(Read_param("Выходные данные: %s",out_file)){
  Read_param("Шаг выходных*: %f",&Dstep);
  OutFile=Err_fopen(out_file,"wt");
 }


 if(Read_param("Входные данные: %s",in_file)){
  if(OutFile &&strcmp(in_file,out_file)==0)serror("Input and output files must be different!\n");
  InFile=Err_fopen(in_file,"rt");
  char str[50];
  fgetline(InFile,str,50);
  float x,y,ymax;
  int Count;
  InType=sscanf(str,"%f %f",&x,&y);
  if(InType<1)serror("Invalid input file format!\n");
  if(InType==1)ymax=x;
  else ymax=y;

  switch(InType){
   case 1:
    Count=1;
    while(!feof(InFile)){
     if(fscanf(InFile,"%f",&y)){
      Count++;
      if(ymax<fabs(y))ymax=fabs(y);
     }
    }
    if(Count<2)serror("Must be at least two data points!\n");
    Dx=1./(Count-1);
   break;
   case 2:
    Count=1;
    Xmin=x;
    while(!feof(InFile)){
     if(fscanf(InFile,"%f %f",&x,&y)){
      Count++;
      if(ymax<fabs(y))ymax=fabs(y);
     }
    }
    XCoeff=1./(x-Xmin);
   break;
  }
  if(ymax<1e-19)serror("The data is too flat!\n");
  YCoeff=1./ymax;
  fseek(InFile,0.,SEEK_SET);
 }

 Close_param_file();


 if(n<1 || n>15)serror("The vavue of p must be between 1 and 15!\n");
 n=2<<(n-1);
 Data=(float *)malloc(2*n*sizeof(float));
 if(!Data)serror("Memory allocation error.\n");
 for(int i=0;i<n;i++)Data[2*i+1]=0;
 CheckStep=1./n;

 int g=EGA,gm=EGAHI;
 initgraph(&g,&gm,"c:\\oldwinap\\borlandc\\bgi");
 pinit(0, 30,30, 600,320, 0,1,-1,1);
 pfield(0, WHITE,BLUE,0,1,1);
 setplotstyle(LIGHTRED,SOLID);
// psetnorepeat(1);
 pdrawplot(FUNCTION,Signal,0.,1.,(CheckStep<0.001 ? CheckStep/4 : 0.00025));
 if(OutFile)fclose(OutFile);
 if(InType)fclose(InFile);

 pfield(0, WHITE,NO_DRAW,0,1,1);
 pgridx(0.5,6);
 pgridx(0.25,3);
 pgridy(1,6);
 pgridy(0.5,3);
 pscalex(0.5,7,"%.1f");
 psetnames("t",UP,"x",RIGHT);
 grprintxy(200,330,"Press any key to continue...");
 getch();
 setfillstyle(SOLID_FILL,0);
 bar(200,330,640,350);

 four(Data,(long)n,1);

 float amax=0,a;
 int imax=0;
 for(i=0;i<n;i++){
  a=sqrt(Data[2*i]*Data[2*i]+Data[2*i+1]*Data[2*i+1]);
  if(a>amax){
   amax=a;
   imax=i;
  }
 }

 ResFreq=(float)imax;
 if(imax!=0){
  ResAcos=(Data[2*imax]+Data[2*n-2*imax])/n;
  ResAsin=(Data[2*imax+1]-Data[2*n-2*imax+1])/n;
 }
 else{
  ResAcos=Data[2*imax]/n;
  ResAsin=Data[2*imax+1]/n;
 }

 setplotstyle(YELLOW,SOLID);
 N=n;
 pdrawplot(FUNCTION,All,0.,1.,0.002);
 setplotstyle(WHITE,SOLID);
 pdrawplot(FUNCTION,ResSignal,0.,1.,0.002);
 setcolor(WHITE);
 if(imax!=0)grprintxy(200,330,"Signal period: %.4f",(float)1./imax);
 else grprintxy(200,330,"Signal period: infinite");
 getch();
 closegraph();
}