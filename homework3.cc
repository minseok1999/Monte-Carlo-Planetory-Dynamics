#include <TFile.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <iostream>
#include <vector>
#include <TVectorD.h>
#include <TVectorT.h>
#include <TAxis.h>
#include <random>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void homework3() 

{

  
    // make a file
    TFile* file = new TFile("homework3.root", "recreate");

   //canvas position and size setting
   TCanvas *c1 = new TCanvas("c1","Radial distance",5000,100,800,600);
   c1->SetGrid();
   
   //number of Iterations
   const Int_t n = 300000;
   
   //Size of sampling points
   const Int_t s=16;

   //defining type
   Double_t r[s], theta[s],tau[s];
  
   //time bin
   const Double_t dt=1.0188*3.14159265/(s-1);

   //Initialization
   //Linear ansatz
   const Double_t rint=0.9;
   const Double_t rfin=1.125;

   for(Int_t i=0; i<s;i++)
   {
     r[i]=rint+i*(rfin-rint)/(s-1);
     theta[i]=3.14159265/(s-1)*i;
     
     tau[i]=i*dt;
     printf(" i %i %f %f \n",i,r[i],theta[i]);
   }

   //defining type of Perturbed array

   Double_t x[s], y[s];
   x[0]=r[0];
   x[s-1]=r[s-1];
   y[0]=theta[0];
   y[s-1]=theta[s-1];
   
   //Size of random fluctuation

   const Double_t rfluc=0.0005;
   const Double_t thetafluc=0.003;

   //Action comparator
   Double_t a[n], b[n];
   
    a[0]=0;
    b[0]=0;
    
    //Monte Carlo Method

   for(Int_t i=0;i<n;i++)
   {
     
     //action for current path calculated as a[i]
     a[i+1]=0;
     for(Int_t j=0;j<s-1;j++)
     {
     a[i+1]= a[i+1]+((pow(((r[j+1]-r[j])/dt),2)+pow(((r[j]+r[j+1])/2*(theta[j+1]-theta[j])/dt),2))/2+(2/(r[j]+r[j+1])))*dt;
     }

     //perturb current path by random number

  
       if(a[i] <= b[i])
       {
           for(Int_t j=1; j<s-1;j++)
           {
        x[j]=  (0.5-(double) rand()/RAND_MAX)* ((double) rand()/RAND_MAX)*rfluc +r[j];
        y[j]=  (0.5-(double) rand()/RAND_MAX)* ((double) rand()/RAND_MAX)*thetafluc +theta[j];
           }
       }

      //Supervise learning --> Step furthered once lowered from last iteration

      if(a[i]>b[i])
       {
        for(Int_t j=1; j<s-1;j++)
           {
         x[j]=x[j]+(x[j]-r[j]);
         y[j]=y[j]+(y[j]-theta[j]);
           }
          }
     

     //action for perturbed path calculated as b[i]
      b[i+1]=0;
      for(Int_t j=0;j<s-1;j++)
     {
     b[i+1]= b[i+1]+((pow(((x[j+1]-x[j])/dt),2)+pow(((x[j]+x[j+1])/2*(y[j+1]-y[j])/dt),2))/2+(2/(x[j]+x[j+1])))*dt; 
     }
   
    //update with lower action in accordance with least action principle

    if(a[i+1]>b[i+1])
    {
      for(Int_t j=0;j<s;j++)
      {
        r[j]=x[j];
        theta[j]=y[j];
      }
    }
     //display current value of action along optimization process 
     printf(" i %i %f %f \n",i,a[i+1],b[i+1]);
   }

  

   //Plotting radial distance
   TGraph *gr = new TGraph(s,tau,r);
   gr->SetLineColor(1);
   gr->SetLineWidth(2);
   gr->SetMarkerColor(1);
   gr->SetMarkerStyle(20);
   gr->SetTitle("Radial distance vs Time");
   gr->GetXaxis()->SetTitle("tau");
   gr->GetYaxis()->SetTitle("zeta");
   gr->Draw("ACP");
  
  
   //Plotting angle
   TGraph *gr2 = new TGraph(s,tau,theta);
   gr2->SetLineColor(1);
   gr2->SetLineWidth(2);
   gr2->SetMarkerColor(1);
   gr2->SetMarkerStyle(20);
   gr2->SetTitle("Angle vs Time");
   gr2->GetXaxis()->SetTitle("tau");
   gr2->GetYaxis()->SetTitle("angle");
   gr2->Draw("ACP");
   

   //Trajectory in Cartesian

    Double_t xx[s], yy[s];
  for(Int_t i=0;i<s;i++)
   {
    xx[i]=r[i]*cos(theta[i]);
    yy[i]=r[i]*sin(theta[i]);
   }

   TGraph *traj = new TGraph(s,xx,yy);
   traj->SetLineColor(1);
   traj->SetLineWidth(2);
   traj->SetMarkerColor(1);
   traj->SetMarkerStyle(20);
   traj->SetTitle("Trajectory in Cartesian");
   traj->GetXaxis()->SetTitle("x coordinate");
   traj->GetYaxis()->SetTitle("y coordinate");
   traj->Draw("ACP");
   // saving datas to file
    gr->Write();
    gr2->Write();
    traj->Write();

    file->Write();

}

  int main()
  {

    homework3();

    return 0;
  }
  
  