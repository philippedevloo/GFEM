#include "TPZFrac2d.h"

TPZFrac2d::TPZFrac2d(TPZVec<REAL> &fracend, TPZVec<REAL> &fracdir, REAL nu, int planestress) {
    fFracend = fracend;
    fFracdir = fracdir;
    fFracAlpha = atan2(fracdir[1], fracdir[0]);
    if(planestress) {
        fKappa = (3-nu)/(1.+nu);
    } else {
        fKappa = 3.-4.*nu;
    }
}

TPZFrac2d::~TPZFrac2d() {
}

#define Cos cos
#define Sin sin
#define Sqrt sqrt
#define kappa fKappa

using namespace std;
void TPZFrac2d::Desloc(const TPZVec<REAL> &x, GFemcolors color, TPZFMatrix<STATE> &phi, TPZFMatrix<STATE> &gradphi) const {
    TPZManVector<REAL,3> xloc(3);
    for(int i=0; i<3; i++) xloc[i] = x[i] - fFracend[i];
    REAL r = sqrt(xloc[0]*xloc[0]+xloc[1]*xloc[1]);
    REAL sqr = sqrt(r);
    REAL theta = atan2(xloc[1],xloc[0]);
    REAL thetafrac = theta - fFracAlpha;
    TPZFNMatrix<4,REAL> rotate(2,2),rotateT(2,2);
    rotate(0,0) = xloc[0]/r;
    rotate(0,1) = xloc[1]/r;
    rotate(1,0) = -xloc[1]/r;
    rotate(1,1) = xloc[0]/r;
    if(IsZero(r)) {
        theta = 0;
        thetafrac = fFracAlpha;
        r = 1.e-10;
        sqr = 1.e-5;
        rotate.Identity();
    }
    rotate.Transpose(&rotateT);
    switch(fSimulationType) {
        case Vector2d:
        {
            TPZFNMatrix<4,REAL> phirt(2,1), gradphirt(2,2), phixy(2,1), gradphixy(2,2);
            REAL u0,u1,gradu0[2],gradu1[2];
            u0 = Sqrt(r)*Cos(thetafrac/2.)*(fKappa - Cos(thetafrac));
            u1 = Sqrt(r)*(-fKappa + Cos(thetafrac))*Sin(thetafrac/2.);
            gradu0[0] = (Cos(thetafrac/2.)*(kappa - Cos(thetafrac)))/(2.*Sqrt(r));
            gradu0[1] = ((2 + kappa + Cos(thetafrac))*Sin(thetafrac/2.))/(2.*Sqrt(r));
            gradu1[0] = ((-kappa + Cos(thetafrac))*Sin(thetafrac/2.))/(2.*Sqrt(r));
            gradu1[1] = (Cos(thetafrac/2.)*(-2 + kappa + Cos(thetafrac)))/(2.*Sqrt(r));
            phirt(0,0) = u0;
            phirt(1,0) = u1;
            gradphirt(0,0) = gradu0[0];
            gradphirt(0,1) = gradu0[1];
            gradphirt(1,0) = gradu1[0];
            gradphirt(1,1) = gradu1[1];
            phixy = rotateT*phirt;
            gradphixy = rotateT*gradphirt*rotate;
            phi(0,0) = phixy(0,0);
            phi(1,0) = phixy(1,0);
            gradphi(0,0) = gradphixy(0,0);
            gradphi(1,0) = gradphixy(0,1);
            gradphi(3,0) = gradphixy(1,0);
            gradphi(4,0) = gradphixy(1,1);
        }
        {
            REAL u0,u1,gradu0[2],gradu1[2];
            TPZFNMatrix<4,STATE> phirt(2,1), gradphirt(2,2), phixy(2,1), gradphixy(2,2);
            u0 = Sqrt(r)*(2 - kappa + 3*Cos(thetafrac))*Sin(thetafrac/2.);
            u1 = -(Sqrt(r)*Cos(thetafrac/2.)*(2 + kappa - 3*Cos(thetafrac)));
            gradu0[0] = ((2 - kappa + 3*Cos(thetafrac))*Sin(thetafrac/2.))/(2.*Sqrt(r));
            gradu0[1] = (Cos(thetafrac/2.)*(kappa + 3*Cos(thetafrac)))/(2.*Sqrt(r));
            gradu1[0] = -0.5*(Cos(thetafrac/2.)*(2 + kappa - 3*Cos(thetafrac)))/Sqrt(r);
            gradu1[1] = -0.5*((kappa + 3*Cos(thetafrac))*Sin(thetafrac/2.))/Sqrt(r);
            phirt(0,0) = u0;
            phirt(1,0) = u1;
            gradphirt(0,0) = gradu0[0];
            gradphirt(0,1) = gradu0[1];
            gradphirt(1,0) = gradu1[0];
            gradphirt(1,1) = gradu1[1];
            phixy = rotateT*phirt;
            gradphixy = rotateT*gradphirt*rotate;
            phi(0,1) = phixy(0,0);
            phi(1,1) = phixy(1,0);
            gradphi(0,1) = gradphixy(0,0);
            gradphi(1,1) = gradphixy(0,1);
            gradphi(3,1) = gradphixy(1,0);
            gradphi(4,1) = gradphixy(1,1);
        }
        break;
    case Scalar2d:
        {
            TPZFNMatrix<4,REAL> phirt(1,1), gradphirt(2,1), gradphixy(2,2);
            REAL u0,gradu0[2];
            u0 = Sin(thetafrac/2.)/Sqrt(r);
            gradu0[0] = -0.5*Sin(theta/2.)/(r*sqr);
            gradu0[1] = Cos(theta/2.)/(2.*r*sqr);
            phirt(0,0) = u0;
            gradphirt(0,0) = gradu0[0];
            gradphirt(1,0) = gradu0[1];
            gradphixy = rotateT*gradphirt;
            phi(0,0) = phirt(0,0);
            gradphi(0,0) = gradphixy(0,0);
            gradphi(1,0) = gradphixy(1,0);
            gradphi(2,0) = 0.;
        }
        break;
    default:
        DebugStop();
    }
}