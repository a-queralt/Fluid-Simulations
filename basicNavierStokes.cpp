#include <fstream>
#include <direct.h>
#include <math.h>
#include <iostream>

using namespace std;

int main(){
    //Declaració variables

    //Generals
    int nR = 40;
    int nC = 160;
    double l=4;
    double h=1;
    double dx;
    double dy;
    double dt=0.0001;
    double convCrit=0.001;
    double itMax=1000;
    double maxTimestep=100000000;
    double rho=7500;
    double steadyState=0.000000001;
    double mu=1;
    double f1=1.5;
    double f2=-0.5;
    bool nordParalel = false;
    bool sudParalel = false;
    bool estParalel = false;
    bool oestParalel = false;
    int saveInterval = 1000;

    //Condicions contorn
    double uN = 1;
    double vN = 0;
    double uS = 0;
    double vS = 0;
    double uE = 0;
    double vE = 0;
    double uW = 0;
    double vW = 0;

    //Malla principal
    double P[nR][nC];
    double posX[nR][nC];
    double posY[nR][nC];
    double velX[nR][nC];
    double velY[nR][nC];    
    double vel[nR][nC];
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
    //càlcul coefficients ai
    for(r=0;r<nR;r++){
        for(c=0;c<nC;c++){
            if(c==0) aw[r][c] = 0;
            else aw[r][c] = 1;
            if(c==nR-1) ae[r][c] = 0;
            else ae[r][c] = 1;
            if(r==0) as[r][c] = 0;
            else as[r][c] = 1;
            if(r==nR-1) an[r][c] = 0;
            else an[r][c] = 1;
            
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
    double un,us,ue,uw, vn,vs,ve,vw;    
    double Raux, auxVel;
    fstream timeSave;

    //Valors Inicials
    for(r=0;r<nR;r++){
        for(c=0;c<nC;c++){
            P[r][c] = 0;
        }
    }

    timeStep=0;
    double timeDiff = 1;
    while(timeStep<maxTimestep /*& steadyState<timeDiff*/){
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

                //Ru i up
                Raux = Ru[r][c];
                Ru[r][c] = -rho*dx*(vn*un-vs*us+ue*ue-uw*uw)+mu*(u[r+1][c]+u[r-1][c]+u[r][c+1]+u[r][c-1]-4*u[r][c]);
                
                up[r][c] = u[r][c]+dt/(rho*pow(dx,2))*(f1*Ru[r][c]+f2*Raux);
                //up[r][c] = u[r][c]+dt/rho*Ru[r][c];
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

                //Rv i vp
                Raux = Rv[r][c];
                Rv[r][c] = -rho*dy*(vn*vn-vs*vs+ue*ve-uw*vw)+mu*(v[r+1][c]+v[r-1][c]+v[r][c+1]+v[r][c-1]-4*v[r][c]);

                vp[r][c] = v[r][c]+dt/(rho*pow(dy,2))*(f1*Rv[r][c]+f2*Raux);
                //vp[r][c] = v[r][c]+dt/rho*Rv[r][c];
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
               /* if(r==0 || r==nR-1 || c==0 || c==nC-1) bp[r][c] = 0;
                else bp[r][c] =  -rho/dt*(vp[r+1][c]*dx-vp[r][c]*dx+up[r][c+1]*dy-up[r][c]*dy);*/
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
                //if(c==0 & r==0) P[r][c] = 0;
                //else{                    
                    oldVal = P[r][c];
                    num = bp[r][c];
                    if(r!=0) num +=P[r-1][c]*as[r][c];
                    if(r!=nR-1) num+= P[r+1][c]*an[r][c];
                    if(c!=0) num+=P[r][c-1]*aw[r][c];
                    if(c!=nC-1) num+=P[r][c+1]*ae[r][c];
                    
                    P[r][c] = num/ap[r][c];

                    currentDiff = fabs(oldVal-P[r][c]);
                    maxDiff=max(maxDiff,currentDiff);
                //}
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
        timeSave<<"_r,_c,posX,posY,posZ,bp,P,velX,velY,vel\n";
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
                    timeSave<<vel[r][c]<<"\n";
                }
            }
        timeSave.close();   
        }
    timeStep++;
    }
}