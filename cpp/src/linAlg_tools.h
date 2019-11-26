#ifndef CRXS__LA_tools_H
#define CRXS__LA_tools_H

#include "math.h"

namespace CRXS {
    
    class LA{
        
        public:
        /// Length of a vector a
        static double len (double* a);
        
        /// cross product of vector v1 and v2 returned into vector cross_res
        static void cross  (double* v1, double* v2, double* cross_res);
        
        /// dot product of vectors v1 and v2
        static double dot(double* v1, double* v2);
        
        /// distance of vector p from the plane through vectors a,b,c
        static double dist( double* a, double *b , double* c, double* p );
        static double sign(double x);
        
        /// true if vector p is inside of the volume of vectors a,b,c,d
        static bool inside(double* a, double *b , double* c, double* d, double* p);
    };
    
    class Integration{
        
        public:
        static double integrate_trapeze( double (*integrand)(double, void*), double min, double max, void* parameter ){
            double res = 0;
            double dd  = (max-min)/steps;
            double d;
            for (int i=0; i<steps; i++) {
                d = min + dd * ( 0.5 + i );
                res += integrand(d,parameter);
            }
            res *= dd;
            return res;
        };
        static int    steps;
        static void   SetTrapezeIntegrationSteps( int _steps ){steps=_steps;};
        
    };
}

#endif




