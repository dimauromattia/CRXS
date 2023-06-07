#include "iostream"
#include "math.h"
#include "xs.h"
#include "xs_definitions.h"

using namespace  CRXS;

int main(){
    
    
    std::cout << "XS_definitions::el_pbarp(10)" << std::endl;
    std::cout <<  XS_definitions::el_pbarp(10)  << std::endl;
    
    std::cout << std::endl;
    std::cout << "Checks on antiproton cross sections" << std::endl;
    std::cout << std::endl;

    std::cout << std::endl;
    std::cout << "Testing inv_AA_pbar_CM for the different CS params " << std::endl;
    //double Tn_proj_LAB, double T_pbar_LAB, double eta_LAB, int A_projectile, int N_projectile, int A_target, int N_target, int parametrization
    std::cout << std::endl;
    std::cout << "XS::inv_AA_pbar_CM(100, 3, 2, 1, 0, 1, 0, 1)" << std::endl;
    std::cout <<  XS::inv_AA_pbar_CM(100, 3, 2, 1, 0, 1, 0, 1)  << std::endl;
    std::cout << "XS::inv_AA_pbar_CM(100, 3, 2, 1, 0, 1, 0, 2)" << std::endl;
    std::cout <<  XS::inv_AA_pbar_CM(100, 3, 2, 1, 0, 1, 0, 2)  << std::endl;
    std::cout << "XS::inv_AA_pbar_CM(100, 3, 2, 1, 0, 1, 0, 3)" << std::endl;
    std::cout <<  XS::inv_AA_pbar_CM(100, 3, 2, 1, 0, 1, 0, 3)  << std::endl;
    std::cout << "XS::inv_AA_pbar_CM(100, 3, 2, 1, 0, 1, 0, 4)" << std::endl;
    std::cout <<  XS::inv_AA_pbar_CM(100, 3, 2, 1, 0, 1, 0, 4)  << std::endl;
    std::cout << "XS::inv_AA_pbar_CM(100, 3, 2, 1, 0, 1, 0, 5)" << std::endl;
    std::cout <<  XS::inv_AA_pbar_CM(100, 3, 2, 1, 0, 1, 0, 5)  << std::endl;
    std::cout << "XS::inv_AA_pbar_CM(100, 3, 2, 1, 0, 1, 0, 10)" << std::endl;
    std::cout <<  XS::inv_AA_pbar_CM(100, 3, 2, 1, 0, 1, 0, 10)  << std::endl;
    std::cout << "XS::inv_AA_pbar_CM(100, 3, 2, 1, 0, 1, 0, 11)" << std::endl;
    std::cout <<  XS::inv_AA_pbar_CM(100, 3, 2, 1, 0, 1, 0, 11)  << std::endl;

    std::cout << std::endl;
    std::cout << "Testing inv_AA_pbar_LAB for the different CS params " << std::endl;
    std::cout << std::endl;
    std::cout << "XS::inv_AA_pbar_LAB(100, 3, 10, 1, 0, 1, 0, 1)" << std::endl;
    std::cout <<  XS::inv_AA_pbar_LAB(100, 3, 10, 1, 0, 1, 0, 1)  << std::endl;
    std::cout << "XS::inv_AA_pbar_LAB(100, 3, 10, 1, 0, 1, 0, 2)" << std::endl;
    std::cout <<  XS::inv_AA_pbar_LAB(100, 3, 10, 1, 0, 1, 0, 2)  << std::endl;
    std::cout << "XS::inv_AA_pbar_LAB(100, 3, 10, 1, 0, 1, 0, 3)" << std::endl;
    std::cout <<  XS::inv_AA_pbar_LAB(100, 3, 10, 1, 0, 1, 0, 3)  << std::endl;
    std::cout << "XS::inv_AA_pbar_LAB(100, 3, 10, 1, 0, 1, 0, 4)" << std::endl;
    std::cout <<  XS::inv_AA_pbar_LAB(100, 3, 10, 1, 0, 1, 0, 4)  << std::endl;
    std::cout << "XS::inv_AA_pbar_LAB(100, 3, 10, 1, 0, 1, 0, 5)" << std::endl;
    std::cout <<  XS::inv_AA_pbar_LAB(100, 3, 10, 1, 0, 1, 0, 5)  << std::endl;
    std::cout << "XS::inv_AA_pbar_LAB(100, 3, 10, 1, 0, 1, 0, 10)" << std::endl;
    std::cout <<  XS::inv_AA_pbar_LAB(100, 3, 10, 1, 0, 1, 0, 10)  << std::endl;
    std::cout << "XS::inv_AA_pbar_LAB(100, 3, 10, 1, 0, 1, 0, 11)" << std::endl;
    std::cout <<  XS::inv_AA_pbar_LAB(100, 3, 10, 1, 0, 1, 0, 11)  << std::endl;

    std::cout << std::endl;
    std::cout << "Testing dE_AA_pbar_LAB for the different CS params " << std::endl;
    std::cout << std::endl;
    //double Tn_proj_LAB, double T_pbar_LAB, int A_projectile, int N_projectile, int A_target, int N_target, int parametrization
    std::cout << "XS::dE_AA_pbar_LAB(100, 5, 1, 0, 1, 0, 1)" << std::endl;
    std::cout <<  XS::dE_AA_pbar_LAB(100, 5, 1, 0, 1, 0, 1)  << std::endl;
    std::cout << "XS::dE_AA_pbar_LAB(100, 5, 1, 0, 1, 0, 2)" << std::endl;
    std::cout <<  XS::dE_AA_pbar_LAB(100, 5, 1, 0, 1, 0, 2)  << std::endl;
    std::cout << "XS::dE_AA_pbar_LAB(100, 5, 1, 0, 1, 0, 3)" << std::endl;
    std::cout <<  XS::dE_AA_pbar_LAB(100, 5, 1, 0, 1, 0, 3)  << std::endl;
    std::cout << "XS::dE_AA_pbar_LAB(100, 5, 1, 0, 1, 0, 4)" << std::endl;
    std::cout <<  XS::dE_AA_pbar_LAB(100, 5, 1, 0, 1, 0, 4)  << std::endl;
    std::cout << "XS::dE_AA_pbar_LAB(100, 5, 1, 0, 1, 0, 5)" << std::endl;
    std::cout <<  XS::dE_AA_pbar_LAB(100, 5, 1, 0, 1, 0, 5)  << std::endl;
    std::cout << "XS::dE_AA_pbar_LAB(100, 5, 1, 0, 1, 0, 10)" << std::endl;
    std::cout <<  XS::dE_AA_pbar_LAB(100, 5, 1, 0, 1, 0, 10)  << std::endl;
    std::cout << "XS::dE_AA_pbar_LAB(100, 5, 1, 0, 1, 0, 11)" << std::endl;
    std::cout <<  XS::dE_AA_pbar_LAB(100, 5, 1, 0, 1, 0, 11)  << std::endl;
    
    std::cout << std::endl;
    //double Tn_proj_LAB, double T_pbar_LAB, int A_projectile, int N_projectile, int A_target, int N_target, int parametrization
    std::cout << "XS::dE_AA_pbar_LAB_incNbarAndHyperon(100, 3)"   << std::endl;
    std::cout <<  XS::dE_AA_pbar_LAB_incNbarAndHyperon(100, 3)    << std::endl;
    std::cout << "XS::dE_AA_p_LAB(100, 3)"   << std::endl;
    std::cout <<  XS::dE_AA_p_LAB(100, 3)    << std::endl;

    std::cout << std::endl;
    std::cout << "Checks on Dbar cross section" << std::endl;
    std::cout << std::endl;
    
    std::cout << std::endl;
    std::cout << "Testing inv_AA_Dbar_LAB for the different antinuclei " << std::endl;
    std::cout << std::endl;
    //double Tn_proj_LAB, double Tn_Hebar_LAB, double eta_LAB, int A_projectile, int N_projectile, int A_target, int N_target, int parametrization, int coalescence
    std::cout << "XS::inv_AA_Dbar_LAB(100, 3, 10, 1, 0, 1, 0, 2, 1, 0.190)" << std::endl;
    std::cout <<  XS::inv_AA_Dbar_LAB(100, 3, 10, 1, 0, 1, 0, 2, 1, 0.190)  << std::endl;
    std::cout << "XS::inv_AA_He3bar_LAB(100, 3, 10, 1, 0, 1, 0, 2, 1, 0.190)" << std::endl;
    std::cout <<  XS::inv_AA_He3bar_LAB(100, 3, 10, 1, 0, 1, 0, 2, 1, 0.190)  << std::endl;
    std::cout << "XS::inv_AA_He4bar_LAB(100, 3, 10, 1, 0, 1, 0, 2, 1, 0.190)" << std::endl;
    std::cout <<  XS::inv_AA_He4bar_LAB(100, 3, 10, 1, 0, 1, 0, 2, 1, 0.190)  << std::endl;

    ////double XS::dEn_AA_He3bar_LAB( double Tn_proj_LAB, double Tn_Hebar_LAB, int A_projectile, int N_projectile, int A_target, int N_target, int parametrization, int coalescence )
    std::cout << std::endl;
    std::cout << "Testing dEn_AA_Dbar_LAB for the two coalescence models " << std::endl;
    std::cout << std::endl;
    std::cout << "XS::dEn_AA_Dbar_LAB(100, 10, 1, 0, 1, 0, 2, 1, 0.190)" << std::endl;
    std::cout <<  XS::dEn_AA_Dbar_LAB(100, 10, 1, 0, 1, 0, 2, 1, 0.190)  << std::endl;
    std::cout << "XS::dEn_AA_Dbar_LAB(100, 10, 1, 0, 1, 0, 2, 2, 0.190)" << std::endl;
    std::cout <<  XS::dEn_AA_Dbar_LAB(100, 10, 1, 0, 1, 0, 2, 2, 0.190)  << std::endl;
    std::cout << "XS::dEn_AA_Dbar_LAB(100, 10, 1, 0, 1, 0, 2, 3, 0.190)" << std::endl;
    std::cout <<  XS::dEn_AA_Dbar_LAB(100, 10, 1, 0, 1, 0, 2, 3, 0.190)  << std::endl;
    std::cout << "XS::dEn_AA_Dbar_LAB(1000, 10, 1, 0, 1, 0, 2, 1, 0.190)" << std::endl;
    std::cout <<  XS::dEn_AA_Dbar_LAB(1000, 10, 1, 0, 1, 0, 2, 1, 0.190)  << std::endl;
    std::cout << "XS::dEn_AA_Dbar_LAB(1000, 10, 1, 0, 1, 0, 2, 2, 0.190)" << std::endl;
    std::cout <<  XS::dEn_AA_Dbar_LAB(1000, 10, 1, 0, 1, 0, 2, 2, 0.190)  << std::endl;
    std::cout << "XS::dEn_AA_Dbar_LAB(1000, 10, 1, 0, 1, 0, 2, 3, 0.190)" << std::endl;
    std::cout <<  XS::dEn_AA_Dbar_LAB(1000, 10, 1, 0, 1, 0, 2, 3, 0.190)  << std::endl;

    std::cout << std::endl;
    std::cout << "Testing dEn_AA_Dbar_LAB for the different coalescence momentum values " << std::endl;
    std::cout << std::endl;
    std::cout << "XS::dEn_AA_Dbar_LAB(100, 10, 1, 0, 1, 0, 2, 1, 0.180)" << std::endl;
    std::cout <<  XS::dEn_AA_Dbar_LAB(100, 10, 1, 0, 1, 0, 2, 1, 0.180)  << std::endl;
    std::cout << "XS::dEn_AA_Dbar_LAB(100, 10, 1, 0, 1, 0, 2, 1, 0.190)" << std::endl;
    std::cout <<  XS::dEn_AA_Dbar_LAB(100, 10, 1, 0, 1, 0, 2, 1, 0.190)  << std::endl;
    std::cout << "XS::dEn_AA_Dbar_LAB(100, 10, 1, 0, 1, 0, 2, 1, 0.200)" << std::endl;
    std::cout <<  XS::dEn_AA_Dbar_LAB(100, 10, 1, 0, 1, 0, 2, 1, 0.200)  << std::endl;
    std::cout << "XS::dEn_AA_Dbar_LAB(100, 10, 1, 0, 1, 0, 2, 1, 0.210)" << std::endl;
    std::cout <<  XS::dEn_AA_Dbar_LAB(100, 10, 1, 0, 1, 0, 2, 1, 0.210)  << std::endl;

    std::cout << std::endl;
    std::cout << "Testing dEn_AA_Dbar_LAB for the different CS parametrizations " << std::endl;
    std::cout << std::endl;
    std::cout << "XS::dEn_AA_Dbar_LAB(100, 10, 1, 0, 1, 0, 1, 1, 0.190)" << std::endl;
    std::cout <<  XS::dEn_AA_Dbar_LAB(100, 10, 1, 0, 1, 0, 1, 1, 0.190)  << std::endl;
    std::cout << "XS::dEn_AA_Dbar_LAB(100, 10, 1, 0, 1, 0, 2, 1, 0.190)" << std::endl;
    std::cout <<  XS::dEn_AA_Dbar_LAB(100, 10, 1, 0, 1, 0, 2, 1, 0.190)  << std::endl;
    std::cout << "XS::dEn_AA_Dbar_LAB(100, 10, 1, 0, 1, 0, 3, 1, 0.190)" << std::endl;
    std::cout <<  XS::dEn_AA_Dbar_LAB(100, 10, 1, 0, 1, 0, 3, 1, 0.190)  << std::endl;
    std::cout << "XS::dEn_AA_Dbar_LAB(100, 10, 1, 0, 1, 0, 4, 1, 0.190)" << std::endl;
    std::cout <<  XS::dEn_AA_Dbar_LAB(100, 10, 1, 0, 1, 0, 4, 1, 0.190)  << std::endl;
    std::cout << "XS::dEn_AA_Dbar_LAB(100, 10, 1, 0, 1, 0, 5, 1, 0.190)" << std::endl;
    std::cout <<  XS::dEn_AA_Dbar_LAB(100, 10, 1, 0, 1, 0, 5, 1, 0.190)  << std::endl;
    std::cout << "XS::dEn_AA_Dbar_LAB(100, 10, 1, 0, 1, 0, 10, 1, 0.190)" << std::endl;
    std::cout <<  XS::dEn_AA_Dbar_LAB(100, 10, 1, 0, 1, 0, 10, 1, 0.190)  << std::endl;
    std::cout << "XS::dEn_AA_Dbar_LAB(100, 10, 1, 0, 1, 0, 11, 1, 0.190)" << std::endl;
    std::cout <<  XS::dEn_AA_Dbar_LAB(100, 10, 1, 0, 1, 0, 11, 1, 0.190)  << std::endl;

    std::cout << "XS::dEn_AA_Dbar_LAB(100, 10, 1, 0, 1, 0, 2, 1, 0.190)" << std::endl;
    std::cout <<  XS::dEn_AA_Dbar_LAB(100, 10, 1, 0, 1, 0, 2, 1, 0.190)  << std::endl;
    std::cout << "XS::dEn_AA_He3bar_LAB(100, 10, 1, 0, 1, 0, 2, 1, 0.190)" << std::endl;
    std::cout <<  XS::dEn_AA_He3bar_LAB(100, 10, 1, 0, 1, 0, 2, 1, 0.190)  << std::endl;
    std::cout << "XS::dEn_AA_He4bar_LAB(100, 10, 1, 0, 1, 0, 2, 1, 0.190)" << std::endl;
    std::cout <<  XS::dEn_AA_He4bar_LAB(100, 10, 1, 0, 1, 0, 2, 1, 0.190)  << std::endl;

    std::cout << "XS::dEn_AA_Dbar_LAB(100, 10, -2, 0, 1, 0, 2, 1, 0.190)" << std::endl;
    std::cout <<  XS::dEn_AA_Dbar_LAB(100, 10, -2, 0, 1, 0, 2, 1, 0.190)  << std::endl;
    std::cout << "XS::dEn_AA_He3bar_LAB(100, 10, -3, 0, 1, 0, 2, 1, 0.190)" << std::endl;
    std::cout <<  XS::dEn_AA_He3bar_LAB(100, 10, -3, 0, 1, 0, 2, 1, 0.190)  << std::endl;
    std::cout << "XS::dEn_AA_He4bar_LAB(100, 10, -4, 0, 1, 0, 2, 1, 0.190)" << std::endl;
    std::cout <<  XS::dEn_AA_He4bar_LAB(100, 10, -4, 0, 1, 0, 2, 1, 0.190)  << std::endl;

    std::cout << "XS::dEn_DbarA_Dbar_LAB(100, 3, 10, 1, 0)" << std::endl;
    std::cout <<  XS::dEn_DbarA_Dbar_LAB(100, 3, 10, 1, 0)  << std::endl;

    std::cout << "XS::p_coal__VonDoetinchen(20)" << std::endl;
    std::cout <<  XS::p_coal__VonDoetinchen(20*20)  << std::endl;
    
    std::cout << "XS::p_coal__VonDoetinchen(90)" << std::endl;
    std::cout <<  XS::p_coal__VonDoetinchen(90*90)  << std::endl;

    std::cout << "XS::p_coal__VonDoetinchen(200)" << std::endl;
    std::cout <<  XS::p_coal__VonDoetinchen(200*200)  << std::endl;

    //for(int i=0; i<1e9; i++){
    //   std::cout <<  XS::dE_AA_p_LAB(  10, 5 ) << std::endl;
    //}

};
