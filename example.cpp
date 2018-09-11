#include "iostream"
#include "math.h"
#include "xs.h"

using namespace  CRXS;

int main(){
    std::cout << "XS::inv_AA_pbar_LAB(100, 3, 10, 1, 0, 1, 0, 2)" << std::endl;
    std::cout <<  XS::inv_AA_pbar_LAB(100, 3, 10, 1, 0, 1, 0, 2)  << std::endl;
    std::cout << "XS::dE_AA_pbar_LAB_incNbarAndHyperon(100, 3)"   << std::endl;
    std::cout <<  XS::dE_AA_pbar_LAB_incNbarAndHyperon(100, 3)    << std::endl;

};
