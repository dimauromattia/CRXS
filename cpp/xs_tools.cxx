#include "xs_tools.h"
#include "xs_definitions.h"

double inv_pp_pbar_CM__Winkler(double s, double E_pbar, double pT_pbar, double* C_array, int len_C_array){
    return CRXS::XS_definitions::inv_pp_pbar_CM__Winkler( s, E_pbar, pT_pbar, &C_array[0], len_C_array );
};
double deltaHyperon(double s, double* C_array, int len_C_array){
    return CRXS::XS_definitions::deltaHyperon( s, &C_array[0], len_C_array );
};

double deltaIsospin(double s,  double* C_array, int len_C_array){
    return CRXS::XS_definitions::deltaIsospin( s, &C_array[0], len_C_array );
};
double inv_pp_pbar_CM__diMauro( double s, double E_pbar, double pT_pbar, double* C_array, int len_C_array ){
    return CRXS::XS_definitions::inv_pp_pbar_CM__diMauro( s, E_pbar, pT_pbar, &C_array[0], len_C_array );
};
double tot_pp__diMauro(double s){
    return CRXS::XS_definitions::tot_pp__diMauro( s );
};
double el_pp__diMauro(double s){
    return CRXS::XS_definitions::el_pp__diMauro( s );
};
double pbar_overlap_function_projectile(double x_F){
    return CRXS::XS_definitions::pbar_overlap_function_projectile( x_F );
};
double pbar_overlap_function_target(double x_F){
    return CRXS::XS_definitions::pbar_overlap_function_target( x_F );
};
