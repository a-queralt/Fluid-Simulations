#include <fstream>
#include <direct.h>
#include <math.h>
#include <iostream>

using namespace std;

double esquema(string metode, double v1, double v2, double v3, double v4, double d12, double d23, double d34){
    if(metode=="CDS") return 0.5*(v2+v3);
    else if(metode=="UDS"){
        if(v2>0) return v2;
        else return v3;
    }
    else if(metode=="SUDS"){
        if(v2>0) return v2*(1+0.5*d23/d12)-v1*0.5*d23/d12;
        else return v3*(1+0.5*d23/d34)-v4*0.5*d23/d34;
    }
    else if(metode=="QUICK"){
        if(v2>0) return -v1*1/4*pow(d23,2)/((d23+d12)*d12)+v2*((d12+0.5*d23)*0.5*d23)/(d12*d23)+v3*((d12+0.5*d23)*(0.5*d23))/((d12+d23)*(d23));
        else return v2*(0.5*d23*(d34+0.5*d23))/(d23*(d23+d34))+v3*(0.5*d23*(d34+0.5*d23))/(d23*d34)-v4*(0.25*d23*d23)/((d23+d34)*d34);
    } 
    else return EXIT_FAILURE;
}

int main(){
    //Declaració variables

    //Geometria
    int nR = 80;
    int nC = 80;
    double l=1;
    double h=1;
    double dx;
    double dy;

    //Numèriques
    double dt=0.0001;
    double dtcx,dtcy,dtd,dtT;
    double convCrit=0.01;
    double itMax=100;
    double maxTimestep=100000000;    
    double steadyState=0.000000001;
    double f1=1.5;
    double f2=-0.5;
    int saveInterval = 1000;

    //Fñisiques
    double rho=100;
    double mu=1;
    double cp = 0.71;
    double lambda = 1;
    double beta = 1;
    double gx = 0;
    double gy = -1000000/(cp*rho*rho);
    double Tref = 0.5;
    double Tini = 0.5;

    //Condicions contorn
    double uN = 0;
    double vN = 0;
    double uS = 0;
    double vS = 0;
    double uE = 0;
    double vE = 0;
    double uW = 0;
    double vW = 0;    
    bool nordParalel = false;
    bool sudParalel = false;
    bool estParalel = false;
    bool oestParalel = false;

    double tempN = 0.5;
    double tempS = 0.5;
    double tempE = 1;
    double tempW = 0;
    bool nordAdiabatic = true;
    bool sudAdiabatic = true;
    bool estAdiabatic = false;
    bool oestAdiabatic = false;

    //Malla principal
    double P[nR][nC];
    double velX[nR][nC];
    double velY[nR][nC];    
    double vel[nR][nC];
    double T[nR][nC];
    double RT[nR][nC];
    double RTold[nR][nC];

    double posX[nR][nC];
    double posY[nR][nC];
    double dE[nR][nC];
    double dW[nR][nC];
    double dN[nR][nC];
    double dS[nR][nC];
    double sN[nR][nC];
    double sS[nR][nC];
    double sW[nR][nC];
    double sE[nR][nC];

    double an[nR][nC];
    double as[nR][nC];
    double ae[nR][nC];
    double aw[nR][nC];
    double ap[nR][nC];
    double bp[nR][nC];

    //Stagg X
    double up[nR][nC+1];
    double u[nR][nC+1];
    double Ru[nR][nC+1];
    
    //Stagg Y
    double vp[nR+1][nC];
    double v[nR+1][nC];
    double Rv[nR+1][nC];   

    //Per bucles
    int r,c,it,timeStep; 

    string metode = "CDS";

    //càlculs geometrics
    dx = l/nC;
    dy = h/nR;
    for(r=0;r<nR;r++){
        for(c=0;c<nC;c++){
            //Posicions
            posX[r][c] = dx*(c+0.5);
            posY[r][c] = dy*(r+0.5); 
        }
    }
    for(r=0;r<nR;r++){
        for(c=0;c<nC;c++){ 
            //Distàncies
            if(r==0) dS[r][c] = 0;
            else dS[r][c] = posY[r][c]-posY[r-1][c];
            if(r==nR-1) dN[r][c] = 0;
            else dN[r][c] = posY[r+1][c]-posY[r][c];
            if(c==0) dW[r][c] = 0;
            else dW[r][c] = posX[r][c]-posX[r][c-1];
            if(c==nC-1) dE[r][c] = 0;
            else dE[r][c] = posX[r][c+1]-posX[r][c];


            //Superfícies
            if(r==0){
                sE[r][c] = dN[r][c];                
                sW[r][c] = dN[r][c];
            }
            else if(r==nR-1){
                sE[r][c] = dS[r][c];
                sW[r][c] = dS[r][c];
            }
            else{
                sE[r][c] = 0.5*(dN[r][c]+dS[r][c]);
                sW[r][c] = 0.5*(dN[r][c]+dS[r][c]);
            }

            if(c==0){
                sN[r][c] = dE[r][c];
                sS[r][c] = dE[r][c];
            }
            else if(c==nC-1){
                sN[r][c] = dW[r][c];
                sS[r][c] = dW[r][c];
            }
            else{
                sN[r][c] = 0.5*(dE[r][c]+dW[r][c]);
                sS[r][c] = 0.5*(dE[r][c]+dW[r][c]);
            }
        }
    }
    //càlcul coefficients ai
    for(r=0;r<nR;r++){
        for(c=0;c<nC;c++){
            if(c==0) aw[r][c] = 0;
            else aw[r][c] = sW[r][c]/dW[r][c];
            if(c==nR-1) ae[r][c] = 0;
            else ae[r][c] = sE[r][c]/dE[r][c];
            if(r==0) as[r][c] = 0;
            else as[r][c] = sS[r][c]/dS[r][c];
            if(r==nR-1) an[r][c] = 0;
            else an[r][c] = sN[r][c]/dN[r][c];
            
            ap[r][c] = an[r][c]+as[r][c]+ae[r][c]+aw[r][c];
        }
    }

    //Impressió geometria+coeficients
    fstream gridV;
    string name = "NS_grid_vertices.txt";
    gridV.open(name, ios::out);         
    gridV<<"_r,_c,posX,posY,posZ,an,as,ae,aw,ap\n";
    for(r = 0; r<nR; r++){
        for(c = 0; c<nC; c++){
                gridV<<r<<",";
                gridV<<c<<",";        
                gridV<<posX[r][c]<<",";
                gridV<<posY[r][c]<<",";
                gridV<<0<<",";
                gridV<<an[r][c]<<",";
                gridV<<as[r][c]<<",";
                gridV<<ae[r][c]<<",";
                gridV<<aw[r][c]<<",";
                gridV<<ap[r][c]<<"\n";
            }
        }
    gridV.close(); 
    
    //Condicions de contorn

    //Malla X
    for(r=0;r<nR;r++){
        u[r][0] = uW;
        u[r][nC] = uE;

        up[r][0] = uW;
        up[r][nC] = uE;
    }
    for(c=0;c<=nC;c++){
        u[0][c] = uS;
        u[nR-1][c] = uN;
        
        up[0][c] = uS;
        up[nR-1][c] = uN;
    }

    //Malla Y
    for(r=0;r<=nR;r++){
        v[r][0] = vW;
        v[r][nC-1] = vE;

        vp[r][0] = vW;
        vp[r][nC-1] = vE;
    }
    for(c=0;c<nC;c++){
        v[0][c] = vS;
        v[nR][c] = vN;
        
        vp[0][c] = vS;
        vp[nR][c] = vN;
    }

    //Simulació
    double un,us,ue,uw, vn,vs,ve,vw, Tn,Ts,Te,Tw;    
    double Raux, auxVel;
    fstream timeSave;
    double sX,sY,a_N,a_S,a_E,a_W;
    //Valors Inicials
    for(r=0;r<nR;r++){
        for(c=0;c<nC;c++){
            P[r][c] = 0;
            //T[r][c] = posX[r][c];
            T[r][c] = Tini;
            RT[r][c] = 0;
            RTold[r][c] = 0;
        }
    }

    //Malla Principal
    for(r=0;r<nR;r++){
        if(!oestAdiabatic) T[r][0] = tempW;
        else T[r][0] = Tini;
        if(!estAdiabatic) T[r][nC-1] = tempE;
        else T[r][nC-1] = Tini;
    }    
    for(c=0;c<nC;c++){
        if(!sudAdiabatic) T[0][c] = tempS;
        else T[0][c] = Tini;
        if(!nordAdiabatic) T[nR-1][c] = tempN;
        else T[nR-1][c] = Tini;
    }

    timeStep=0;
    double timeDiff = 1;
    double dtAux;
    while(timeStep<maxTimestep /*& steadyState<timeDiff*/){
        //Malla Principal

        for(r=1;r<nR-1;r++){
            for(c=1;c<nC-1;c++){
                Tn = 0.5*(T[r+1][c]+T[r][c]);
                Ts = 0.5*(T[r-1][c]+T[r][c]);
                Te = 0.5*(T[r][c+1]+T[r][c]);
                Tw = 0.5*(T[r][c-1]+T[r][c]);

                RTold[r][c] = RT[r][c];
                RT[r][c] = -rho*cp*(u[r][c+1]*sE[r][c]*Te-u[r][c]*sW[r][c]*Tw+v[r+1][c]*sN[r][c]*Tn-v[r][c]*sS[r][c]*Ts)+lambda*(an[r][c]*T[r+1][c]+as[r][c]*T[r-1][c]+ae[r][c]*T[r][c+1]+aw[r][c]*T[r][c-1]-ap[r][c]*T[r][c]);
            }
        }
        for(r=1;r<nR-1;r++){
            for(c=1;c<nC-1;c++){
                T[r][c] = T[r][c] + dt/(rho*cp*dx*dy)*(f1*RT[r][c]+f2*RTold[r][c]);   
            }
        }
        for(r=0;r<nR;r++){
            if(oestAdiabatic) T[r][0] = T[r][1];
            if(estAdiabatic) T[r][nC-1] = T[r][nC-2];
        }    
        for(c=0;c<nC;c++){
            if(sudAdiabatic) T[0][c] = T[1][c];
            if(nordAdiabatic) T[nR-1][c] = T[nR-2][c];
        }
        
        
        //Malla x
        for(r=1;r<nR-1;r++){
            for(c=1;c<nC;c++){
                //Càlcul de velocitats cares node X
                un = 0.5*(u[r+1][c]+u[r][c]);
                us = 0.5*(u[r-1][c]+u[r][c]);
                ue = 0.5*(u[r][c+1]+u[r][c]);
                uw = 0.5*(u[r][c-1]+u[r][c]);

                vn = 0.5*(v[r+1][c-1]+v[r+1][c]);
                vs = 0.5*(v[r][c-1]+v[r][c]);

                sX = sW[r][c]; 
                sY = dW[r][c];

                if(dN[r][c]!=0) a_N = sY/dN[r][c];
                if(dS[r][c]!=0) a_S = sY/dS[r][c];
                if(dE[r][c]!=0) a_E = sX/dE[r][c];
                if(dW[r][c]!=0) a_W = sX/dW[r][c];

                a_N = sY/dN[r][c];
                a_S = sY/dS[r][c];
                a_E = sX/dE[r][c];
                a_W = sX/dW[r][c];

                //Ru i up
                Raux = Ru[r][c];
                Ru[r][c] = -rho*dx*(vn*un-vs*us+ue*ue-uw*uw)+mu*(u[r+1][c]+u[r-1][c]+u[r][c+1]+u[r][c-1]-4*u[r][c])+gx*rho*beta*(T[r][c]-Tref)*dN[r][c]*dE[r][c];
                

                up[r][c] = u[r][c]+dt/(rho*pow(dx,2))*(f1*Ru[r][c]+f2*Raux);
            }
            if(estParalel){                
                up[r][0] = up[r][1];
            }
            if(oestParalel){
                up[r][nC] = up[r][nC-1];
            }
        }

        //Malla Y        
        for(c=1;c<nC-1;c++){
            for(r=1;r<nR;r++){
                vn = 0.5*(v[r+1][c]+v[r][c]);
                vs = 0.5*(v[r-1][c]+v[r][c]);
                ve = 0.5*(v[r][c+1]+v[r][c]);
                vw = 0.5*(v[r][c-1]+v[r][c]);

                ue = 0.5*(u[r][c+1]+u[r-1][c+1]);
                uw = 0.5*(u[r][c]+u[r-1][c]);

                sX = dS[r][c]; 
                sY = 0.5*(dE[r][c]+dW[r][c]);

                a_N = sY/dN[r][c];
                a_S = sY/dS[r][c];
                a_E = sX/dE[r][c];
                a_W = sX/dW[r][c];

                //Rv i vp
                Raux = Rv[r][c];
                Rv[r][c] = -rho*dy*(vn*vn-vs*vs+ue*ve-uw*vw)+mu*(v[r+1][c]+v[r-1][c]+v[r][c+1]+v[r][c-1]-4*v[r][c])+gy*rho*beta*(T[r][c]-Tref)*dN[r][c]*dE[r][c];

                vp[r][c] = v[r][c]+dt/(rho*pow(dy,2))*(f1*Rv[r][c]+f2*Raux);
            }
            if(nordParalel){
                vp[nR][c] = vp[nR-1][c];
            }
            if(sudParalel){
                vp[0][c] = vp[1][c];
            }
        }

        //Malla Principal
        for(r=0;r<nR;r++){
            for(c=0;c<nC;c++){
                bp[r][c] =  -rho/dt*(vp[r+1][c]*dx-vp[r][c]*dx+up[r][c+1]*dy-up[r][c]*dy);
            }
        }

    //Mètode iteratiu
    double oldVal;
    double currentDiff;
    double maxDiff = 1;
    it = 0;
    double num;
    while(it<itMax & maxDiff>convCrit){
        maxDiff=0;
        for(r=0;r<nR;r++){
            for(c=0;c<nC;c++){             
                    oldVal = P[r][c];
                    num = bp[r][c];
                    if(r!=0) num +=P[r-1][c]*as[r][c];
                    if(r!=nR-1) num+= P[r+1][c]*an[r][c];
                    if(c!=0) num+=P[r][c-1]*aw[r][c];
                    if(c!=nC-1) num+=P[r][c+1]*ae[r][c];
                    
                    P[r][c] = num/ap[r][c];

                    currentDiff = fabs(oldVal-P[r][c]);
                    maxDiff=max(maxDiff,currentDiff);                
            }
        }
        if(it%10000==0 & it!=0) cout<<"it: "<<it<<"\n"<<"Diff max:"<<maxDiff<<"\n\n";
        it++;
    }

    //Noves velocitats

    //Malla X
    timeDiff = 0;
    for(r=1;r<nR-1;r++){
        for(c=1;c<nC;c++){
            auxVel = u[r][c];
            u[r][c] = up[r][c]-dt/rho*(P[r][c]-P[r][c-1])/dx;            
            timeDiff = max(timeDiff, fabs(u[r][c]-auxVel));  
        }
    }
    if(nordParalel){
        for(c=0;c<=nC;c++){
            u[nR-1][c] = u[nR-2][c];
            up[nR-1][c] = u[nR-2][c];
        }
    }
    if(sudParalel){
        for(c=0;c<=nC;c++){
            u[1][c] = u[0][c];
            up[1][c] = u[0][c];
        }
    }    
    if(estParalel){
        for(r=0;r<nR;r++){
            u[r][nC] = u[r][nC-1];
            up[r][nC] = u[r][nC-1];
        }
    }    
    if(oestParalel){
        for(r=0;r<nR;r++){
            u[r][0] = u[r][1];            
            up[r][0] = u[r][1];
        }
    } 
    
    //Malla Y
    for(r=1;r<nR;r++){
        for(c=1;c<nC-1;c++){
            auxVel = v[r][c];
            v[r][c] = vp[r][c]-dt/rho*(P[r][c]-P[r-1][c])/dy;
            timeDiff = max(timeDiff, fabs(v[r][c]-auxVel));  
        }
    }
    if(nordParalel){
        for(c=0;c<nC;c++){
            v[nR][c] = v[nR-1][c];
            vp[nR][c] = v[nR-1][c];
        }
    }
    if(sudParalel){
        for(c=0;c<nC;c++){
            v[1][c] = v[0][c];
            vp[1][c] = v[0][c];
        }
    }    
    if(estParalel){
        for(r=0;r<=nR;r++){
            v[r][nC-1] = v[r][nC-2];
            vp[r][nC-1] = v[r][nC-2];
        }
    }    
    if(oestParalel){
        for(r=0;r<=nR;r++){
            v[r][0] = v[r][1];
            vp[r][0] = v[r][1];
        }
    }
    if(timeStep%saveInterval==0){     
        cout<<"Time diff:"<<timeDiff<<"\n";
        //Velocitats a la malla
        for(r=0;r<nR;r++){
            for(c=0;c<nC;c++){
                velX[r][c] = 0.5*(u[r][c]+u[r][c+1]);
                velY[r][c] = 0.5*(v[r][c]+v[r+1][c]);
                vel[r][c] = sqrt(pow(velX[r][c],2)+pow(velY[r][c],2));
            }
        }

        //Desa timestepg 
        name = "results"+to_string(timeStep)+".txt";
        timeSave.open(name, ios::out);         
        timeSave<<"_r,_c,posX,posY,posZ,bp,P,velX,velY,vel,RT,T\n";
        for(r = 0; r<nR; r++){
            for(c = 0; c<nC; c++){
                    timeSave<<r<<",";
                    timeSave<<c<<",";        
                    timeSave<<posX[r][c]<<",";
                    timeSave<<posY[r][c]<<",";
                    timeSave<<0<<",";
                    timeSave<<bp[r][c]<<",";
                    timeSave<<P[r][c]<<",";
                    timeSave<<velX[r][c]<<",";
                    timeSave<<velY[r][c]<<",";
                    timeSave<<vel[r][c]<<",";
                    timeSave<<RT[r][c]<<",";
                    timeSave<<T[r][c]<<"\n";
                }
            }
        timeSave.close();   
        }
    timeStep++;
    }
    
}