#ifndef TPZFRACTURE2D_H
#define TPZFRACTURE2D_H

#include "pzmanvector.h"
#include "pzreal.h"
#include "pzvec.h"
#include "pzfmatrix.h"

enum GFemcolors {Nocolor, Black, White, Grey};

enum GFemsimulation {Nosimulation, Scalar2d, Vector2d};

class TPZFrac2d {
public:
    TPZFrac2d() {}
    TPZFrac2d(TPZVec<REAL> &fracend, TPZVec<REAL> &fracdir, REAL nu, int planestress = 1);
    ~TPZFrac2d();

    TPZFrac2d(const TPZFrac2d &copy) : fFracend(copy.fFracend), fFracdir(copy.fFracdir), fFracAlpha(copy.fFracAlpha), fKappa(copy.fKappa), fSimulationType(copy.fSimulationType) {}

    TPZFrac2d &operator=(const TPZFrac2d &copy)
    {
        fFracend = copy.fFracend;
        fFracdir = copy.fFracdir;
        fFracAlpha = copy.fFracAlpha;
        fKappa = copy.fKappa;
        fSimulationType = copy.fSimulationType;
        return *this;
    }

    void SetFracData(TPZVec<REAL> &fracend, TPZVec<REAL> &fracdir, REAL nu, GFemsimulation simulation, int planestress = 1)
    {
        fFracend = fracend;
        fFracdir = fracdir;
        fFracAlpha = atan2(fracdir[1], fracdir[0]);
        if(planestress) {
            fKappa = (3-nu)/(1.+nu);
        } else {
            fKappa = 3.-4.*nu;
        }
        fSimulationType = simulation;
    }

    void Desloc(const TPZVec<REAL> &x, GFemcolors color, TPZFMatrix<STATE> &phi, TPZFMatrix<STATE> &gradphi) const;

    std::function<void (const TPZVec<REAL> &x, GFemcolors color, TPZFMatrix<STATE> &phi, TPZFMatrix<STATE> &gradphi)> DeslocFunction() const
    {
        return [this](const TPZVec<REAL> &x, GFemcolors color, TPZFMatrix<STATE> &phi, TPZFMatrix<STATE> &gradphi)
        {
            this->Desloc(x,color,phi,gradphi);
        };
    }

private:
    TPZManVector<REAL,3> fFracend;

    TPZManVector<REAL,3> fFracdir;

    /// @brief Angle of the fracture computed in functions of fFracdir
    REAL fFracAlpha;

    /// @brief Coefficient of the fracture displacement as a function of the poisson coefficient and stress state
    REAL fKappa;

    // Type of simulation : either scalar (Darcy2d) or vector (Elastic 2d)
    GFemsimulation fSimulationType;

};
#endif