#include "iostream"
#include "math.h"
#include "xs.h"

#include "gsl_integration.h"

using namespace  CRXS;

double phi(double x){
    
    if (x<5) {
        return pow(x, -1.0);
    }
    return pow(x, -2.8);
    
}


double integrant(double log_Tn, void* parameters ){
    double Tn = exp(log_Tn);
    double Tpbar = * (double * ) parameters;
    return  Tn*XS::dE_AA_pbar_LAB_incNbarAndHyperon(Tn,Tpbar)*phi(Tn);
}




int main(){
    std::cout << XS::inv_AA_pbar_LAB(100, 3, 10, 1, 0, 1, 0, 2) << std::endl;
    std::cout << XS::dE_AA_pbar_LAB_incNbarAndHyperon(100,3) << std::endl;

//    for (double dTpbar = -1; dTpbar<=4; dTpbar+=1./30) {
//        for (double dT =  0; dT<=6; dT+=1./30) {
//            double Tn = pow(10,dT);
//            double Tpbar = pow(10,dTpbar);
//            double res = XS::dE_AA_pbar_LAB_incNbarAndHyperon(Tn,Tpbar);
//        }
//    }
    
    gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);
    gsl_function F;
    F.function = &integrant;
    double epsabs = 0;
    double epsrel = 1e-4;
    double res, err;
    
    for (double dTpbar = -1; dTpbar<=4; dTpbar+=1./30) {
        double Tpbar = pow(10,dTpbar);
        F.params = &Tpbar;
        gsl_integration_qag(&F, 0, 30, epsabs, epsrel, 1000, 2, w, &res, &err);
        std::cout << res << std::endl;
    }
    
    
    
    
    
    
};
