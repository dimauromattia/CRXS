#include "linAlg_tools.h"

namespace CRXS {
    
    
    
    /// Length of a vector a
    double LA::len (double* a){
        return sqrt(a[0]*a[0]+a[1]*a[1]+a[2]*a[2]);
    }
    
    /// cross product of vector v1 and v2 returned into vector cross_res
    void LA::cross  (double* v1, double* v2, double* cross_res){
        cross_res[0] = v1[1]*v2[2]-v1[2]*v2[1];
        cross_res[1] = v1[2]*v2[0]-v1[0]*v2[2];
        cross_res[2] = v1[0]*v2[1]-v1[1]*v2[0];
    }
    
    /// dot product of vectors v1 and v2
    double LA::dot(double* v1, double* v2){
        return v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2];
    }
    
    /// distance of vector p from the plane through vectors a,b,c
    double LA::dist( double* a, double *b , double* c, double* p ){
        
        double v1 [3] = { a[0]-b[0], a[1]-b[1], a[2]-b[2] };
        double v2 [3] = { a[0]-c[0], a[1]-c[1], a[2]-c[2] };
        double v3 [3] = { a[0]-p[0], a[1]-p[1], a[2]-p[2] };
        
        double n  [3] = {0,0,0};
        cross(v1, v2, n);
        double dn = len(n);
        n[0]      = n[0]/dn;
        n[1]      = n[1]/dn;
        n[2]      = n[2]/dn;
        
        return dot(n, v3);
        
    }
    double LA::sign(double x){
        if (x<0)
        return -1;
        return 1;
    }
    
    /// true if vector p is inside of the volume of vectors a,b,c,d
    bool LA::inside(double* a, double *b , double* c, double* d, double* p){
        double dd, dp;
        
        dd = dist(a, b, c, d);
        if ( fabs(dd)<1e-10 || dd!=dd){
            return false;
        }
        dp = dist(a, b, c, p);
        if ( fabs(dp)<1e-10 || dp!=dp){
            return false;
        }
        dp = dp * sign(dd);
        dd = dd * sign(dd);
        if (dp<0 || dp>dd){
            return false;
        }
        
        dd = dist(a, b, d, c);
        if ( fabs(dd)<1e-10 || dd!=dd){
            return false;
        }
        dp = dist(a, b, d, p);
        if ( fabs(dp)<1e-10 || dp!=dp){
            return false;
        }
        dp = dp * sign(dd);
        dd = dd * sign(dd);
        if (dp<0 || dp>dd){
            return false;
        }
        
        dd = dist(a, c, d, b);
        if ( fabs(dd)<1e-10 || dd!=dd){
            return false;
        }
        dp = dist(a, c, d, p);
        if ( fabs(dp)<1e-10 || dp!=dp){
            return false;
        }
        dp = dp * sign(dd);
        dd = dd * sign(dd);
        if (dp<0 || dp>dd){
            return false;
        }
        
        dd = dist(b, c, d, a);
        if ( fabs(dd)<1e-10 || dd!=dd){
            return false;
        }
        dp = dist(b, c, d, p);
        if ( fabs(dp)<1e-10 || dp!=dp){
            return false;
        }
        dp = dp * sign(dd);
        dd = dd * sign(dd);
        if (dp<0 || dp>dd){
            return false;
        }
        
        
        return true;
    }
    
    int Integration::steps = 1000;
    
}



