#define _USE_MATH_DEFINES
#include <math.h>
#include <stdio.h>
#include "fftw.h"
#include "TVirtualFFT.h"
//#include <complex.h>
#include <complex>

using namespace std; 
typedef complex<double> com;
const com i(0.0, 1.0);

int N = 4*1024;
double Pi = M_PI;
double a = -60e-4;
double b =  60e-4;
double L = 10.0;
double dx = 1.0 * (b - a) / N;

        ///////und parameters///////
double K = 1;
double lamda =10e-9; 
double e = 1;
double c = 3.0e8;

double w0 = 2. * Pi * c / lamda;
double g = 1; 
double L0 = 1;
double k0 = (2.0 * Pi) / lamda;
double sigma = 20e-6;
double W0 = sigma;
        /////helper function//////

void space(double x[N])
{
    for(int n = 0; n < N; n++)
    {
        x[n] = a + n * dx; 
    }
}

void inverce_space(double k[N])
{
    for(int n = 0; n < N/2; n++)
    {
        k[n] = 2 * Pi * n / (b - a);
        k[n + N/2] = 2 * Pi * (n - N/2) / (b - a);
    }
}
        ///////distributions////////

com dfl_dist_und(double x)
{
    return (i*((K*w0*e)/(2*c*c*g))* (Pi - 2)); 
}

void gausian_beam(com dfl[N], double x0)
{ 
    double x[N];
    for(int n = 0; n < N; n++)
    {
        x[n] = a + n * dx; 
        dfl[n] = (exp(-(x[n] - x0)*(x[n] - x0)/(W0*W0)));
    }
}

com H(double z, double k)
{
    std::complex<double > arg = (1.0 - (k * k) / (k0 * k0));
    //return ( exp(i * k0 * z) * exp(-(i * z * k * k)/(2 * k0)));
    return ( exp(i * k0 * z * std::pow(arg, 0.5)));
}
com T(double f, double x)
{
    return (exp(-(i * k0)*(x*x)/(2. * f)));
}
com parity (double x, double y)
{
    return (x + i*y);
}

        /////////DFT///////////

void fftw(com Ft[N], com Fw[N], string D)
{   
    double ReFt[N]; 
    double ImFt[N];
    std::string dir;
    for(int n = 0; n < N; n++)
    {
        ReFt[n] = (double)(real(Ft[n]));
        ImFt[n] = (double)(imag(Ft[n]));
    }
    if(D == "FORWARD")
        dir = "C2CFORWARD";
    if(D == "BACKWARD")
        dir = "C2CFORWARD";
    
    char * writable = new char[dir.size() + 1];
    std::copy(dir.begin(), dir.end(), writable);
    writable[dir.size()] = '\0'; // don't forget the terminating 0

    // don't forget to free the string after finished using it
    TVirtualFFT *fftr2c = TVirtualFFT::FFT(1, &N, writable);//C2CFORWARD
    
    fftr2c->SetPointsComplex(ReFt, ImFt);
    fftr2c->Transform();
    double re, im;
    
    if(D == "BACKWARD")
    {
        for(int n = 0; n < N; n++)
        {
            fftr2c->GetPointComplex(n, re, im);  
            re = (double)(re/N); 
            im = (double)(im/N);
            Fw[n] = parity(re, im);
        }
    }
    if(D == "FORWARD")
    {
        for(int n = 0; n < N; n++)
        {
            fftr2c->GetPointComplex(n, re, im);  
            re = (double)(re); 
            im = (double)(im);
            Fw[n] = parity(re, im);
        }
    }
    delete[] writable;
}

void separate(complex<double>* X, int Num)
{
    complex<double>* b = new complex<double>[Num/2];
    for(int i = 0; i < Num/2; i++)
        b[i] = X[i*2 + 1];
    for(int i = 0; i < Num/2; i++)
        X[i] = X[i*2];
    for(int i = 0; i < Num/2; i++)
        X[i+Num/2] = b[i];
    delete[] b;
}


void My_FFT(com* Y, com* X, int Num, string D)
{  
    for(int i = 0; i < N; i++)
        X[i] = Y[i];

    com w; 
    if(Num < 2){}
    else
    {
        separate(X, Num);
        My_FFT(X, X, Num/2, D);
        My_FFT(X + Num/2, X + Num/2, Num/2, D);
        for(int k = 0; k < Num/2; k++)
        {
            com e = X[k      ]; 
            com o = X[k + Num/2];
           
            if(D == "BACKWARD")
                w = exp(com(0, 2 * Pi * k / Num));
            if(D == "FORWARD")
                w = exp(com(0, -2 * Pi * k / Num));
            
            X[k]          = e + w * o;
            X[k + Num/2]  = e - w * o;
        }
    }
    
}

void FFT_2(com* Y, com* X, int Num, string D)
{
    My_FFT(Y, X, Num, D);
    if(D == "BACKWARD")
    {
        for(int i = 0; i < N; i++)
        {
            X[i] = (com(1./N, 0))*X[i];
        }
    }
}

            ////////draw////////
void Spec(com Fw[N], double Sp[N])
{
    for(int n = 0; n < N; n++)
    {
        Sp[n] = sqrt(real(Fw[n]) * real(Fw[n]) + imag(Fw[n]) * imag(Fw[n]));
        //Sp[n] = real(Fw[n]);
        //Sp[n] = imag(Fw[n]);
    }
}
            //////optics//////
void prop(double z, com dfl[N], string FT)
{
    com FT_dfl[N];
    double x[N], k[N];
    
    if(FT == "MY")
        FFT_2(dfl, FT_dfl, N, "FORWARD");
    if(FT == "FFTW")
        fftw(dfl, FT_dfl, "FORWARD");
 
    inverce_space(k); 
    
    for(int n = 0; n < N; n++)
    {   
        FT_dfl[n] = FT_dfl[n] * H(z, k[n]);
        
    } 
    if(FT == "MY")
        FFT_2(FT_dfl, dfl, N, "BACKWARD");
    if(FT == "FFTW")
        fftw(FT_dfl, dfl, "BACKWARD");
}

void lens(double f, com dfl[N])
{
    double x[N];
    space(x);     
    for(int n = 0; n < N; n++)
    {   
        dfl[n] = dfl[n] * T(f, x[n]);    
    }
}
void sum_dfl(com* dfl, com* dfl1, com* dfl2)
{
    for(int n = 0; n < N; n++)
    {   
       dfl[n] = dfl1[n] + dfl2[n]; 
    }
}

    ////////////main//////////

void projectV1_1()
{ 
        ////// fields ////// 
    //cout << c << endl;
    double x[N], k[N], Sp[N]; 
    com FTGaus[N];
    com Gaus[N];
    com dfl1[N], dfl2[N], dfl[N]; com FTdfl[N]; 
      
        ///// domains ///////////
    space(x);
    inverce_space(k);
    
    int Num = N;
        ////// optic system ///// 
    
    double l = 14;
    double f = 7;
    double x0 = 10e-5; 
    gausian_beam(dfl1, -x0);
    gausian_beam(dfl2,  x0);

    Spec(dfl1, Sp);
    TCanvas* c3 = new TCanvas("c3", "Two sources", 500, 500);
    TGraph* gr4 = new TGraph(N, x, Sp);
    gr4->SetName("without Hun`s");
    gr4->Draw("");
    gr4->SetMarkerStyle(8); 
    gr4->SetMarkerSize(0.5);
    gr4->SetMarkerColor(4);

    Spec(dfl2, Sp);
    //TCanvas* c1 = new TCanvas("c1", "without Hun`s window", 500, 500);
    TGraph* gr5 = new TGraph(N, x, Sp);
    gr5->SetName("without Hun`s window");
    gr5->Draw("pl same");
    gr5->SetMarkerStyle(8); 
    gr5->SetMarkerSize(0.5);
    gr5->SetMarkerColor(6);            
    
    prop(l, dfl1, "MY");
    prop(l, dfl2, "FFTW");
    lens(f, dfl1);
    lens(f, dfl2);
    
    double d = f + f*f*(l - f)/((l - f)*(l - f) + (k0*W0*W0/2)*(k0*W0*W0/2));
    cout << d << endl;
    
    prop(d, dfl1, "MY");
    prop(d, dfl2, "FFTW");

    //prop(10, dfl1, "MY");
    //prop(10, dfl2, "FFTW");

    //sum_dfl(dfl, dfl1, dfl2);
    
    
    ///////DRAW///////
/*
    Spec(dfl2, Sp);
    //TCanvas* c1 = new TCanvas("c1", "without Hun`s window", 500, 500);
    TGraph* gr1 = new TGraph(N, x, Sp);
    gr1->SetName("without Hun`s window");
    gr1->Draw("pl same");
    gr1->SetMarkerStyle(8); 
    gr1->SetMarkerSize(0.5);
    gr1->SetMarkerColor(6);            
    */
    Spec(dfl1, Sp);
    TCanvas* c1 = new TCanvas("c1", "interference partern", 500, 500);
    TGraph* gr3 = new TGraph(N, x, Sp);
    gr3->SetName("without Hun`s window");
    gr3->Draw("");
    gr3->SetMarkerStyle(8); 
    gr3->SetMarkerSize(0.5);
    gr3->SetMarkerColor(4);            

}

