//
// Created by Gustavo Batistela on 5/13/21.
//

#include "TPZGFemDarcyFlow.h"
#include "pzaxestools.h"
// #include "TPZLapack.h"

TPZGFemDarcyFlow::TPZGFemDarcyFlow() : TPZRegisterClassId(&TPZGFemDarcyFlow::ClassId),
                               TPZMatCombinedSpacesT<STATE>(), TPZDarcyFlow() {}

TPZGFemDarcyFlow::TPZGFemDarcyFlow(int id, int dim) : TPZRegisterClassId(&TPZGFemDarcyFlow::ClassId),
                                        TPZDarcyFlow(id,dim) {}



/*
int TPZGFemDarcyFlow::VariableIndex(const std::string &name) const {

    if (!strcmp("Solution", name.c_str())) return 1;
    if (!strcmp("Pressure", name.c_str())) return 1;
    if (!strcmp("Derivative", name.c_str())) return 2;
    if (!strcmp("GradU", name.c_str())) return 2;
    if (!strcmp("KDuDx", name.c_str())) return 3;
    if (!strcmp("KDuDy", name.c_str())) return 4;
    if (!strcmp("KDuDz", name.c_str())) return 5;
    if (!strcmp("NormKDu", name.c_str())) return 6;
    if (!strcmp("MinusKGradU", name.c_str())) return 7;
    if (!strcmp("Flux", name.c_str())) return 7;
    if (!strcmp("POrder", name.c_str())) return 8;
    if (!strcmp("ExactPressure", name.c_str())) return 9;
    if (!strcmp("ExactSolution", name.c_str())) return 9;
    if (!strcmp("ExactFlux", name.c_str())) return 10;
    if (!strcmp("Div", name.c_str())) return 11;
    if (!strcmp("Divergence", name.c_str())) return 11;
    if (!strcmp("ExactDiv", name.c_str())) return 12;
    if (!strcmp("ExactDivergence", name.c_str())) return 12;
    if (!strcmp("FluxL2", name.c_str())) return 13;

    return TPZDarcyFlow::VariableIndex(name);
}

int TPZGFemDarcyFlow::NSolutionVariables(int var) const {

    if (var == 1) return 1;      // Solution/Pressure
    if (var == 2) return fDim;   // Derivative/GradU
    if (var == 3) return 1;      // KDuDx;
    if (var == 4) return 1;      // KDuDy;
    if (var == 5) return 1;      // KDuDz;
    if (var == 6) return 1;      // NormKDu;
    if (var == 7) return fDim;   // MinusKGradU/Flux;
    if (var == 8) return 1;      // POrder
    if (var == 9) return 1;      // ExactPressure/ExactSolution
    if (var == 10) return fDim;  // ExactFlux
    if (var == 11) return 1;     // Div/Divergence
    if (var == 12) return 1;     // ExactDiv/ExactDivergence
    if (var == 13) return fDim;  // FluxL2

    return TPZDarcyFlow::NSolutionVariables(var);
}
*/


int TPZGFemDarcyFlow::ClassId() const {
    return Hash("TPZGFemDarcyFlow") ^ TPZDarcyFlow::ClassId() << 1;
}

TPZMaterial *TPZGFemDarcyFlow::NewMaterial() const {
    return new TPZGFemDarcyFlow(*this);
}

void TPZGFemDarcyFlow::Print(std::ostream &out) const {
    out << "Material Name: " << this->Name() << "\n";
    out << "Material Id: " << TPZDarcyFlow::Id() << "\n";
    out << "Dimension: " << TPZDarcyFlow::Dimension() << "\n\n";
}

/** @name Contribute */
/** @{ */
/**
 * @brief It computes a contribution to the stiffness matrix
 * and load vector at one integration point.
 * @param[in] datavec stores all input data
 * @param[in] weight is the weight of the integration rule
 * @param[out] ek is the element matrix
 * @param[out] ef is the rhs vector
 */
void TPZGFemDarcyFlow::Contribute(const TPZVec<TPZMaterialDataT<STATE>> &datavec,
                        REAL weight,TPZFMatrix<STATE> &ek,
                        TPZFMatrix<STATE> &ef)
{
    TPZVec<REAL>  &x = datavec[1].x;

    STATE fXfLoc = 0;

    if(fForcingFunction) {            // phi(in, 0) = phi_in
        TPZManVector<STATE,1> res(1);
        fForcingFunction(x, res);
        fXfLoc = res[0];
    }

    STATE KPerm = GetPermeability(datavec[0].x);
    int64_t nvec = datavec.size();
    int64_t firsti = 0;
    for(int i = 0; i<nvec; i++) {
        int nshapei = datavec[i].phi.Rows();
        if(nshapei < 1) continue;
        if(datavec[i].fShapeType != TPZMaterialData::EScalarShape) DebugStop();
        int64_t firstj = 0;
        for(int j = 0; j<nvec; j++) {
            int nshapej = datavec[j].phi.Rows();
            if(nshapej < 1) continue;
            ek.AddContribution(firsti,firstj,datavec[i].dphix,true,datavec[j].dphix,false,weight*KPerm);
            firstj += datavec[j].phi.Rows();
        }
        auto nshape = datavec[i].phi.Rows();
        for(int k = 0; k<nshape; k++) ef(firsti+k,0) += (weight*fXfLoc)*datavec[i].phi(k,0);
        firsti += datavec[i].phi.Rows();
    }

}

void TPZGFemDarcyFlow::ContributeBC(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCondT<STATE> &bc)
{
    auto &x = datavec[0].x;
    STATE v2;
    v2 = bc.Val2()[0];

    if(bc.HasForcingFunctionBC()) {            // phi(in, 0) = phi_in                          // JORGE 2013 01 26
        TPZManVector<STATE> res(1);
        TPZFNMatrix<3,STATE> dres(3,1);
        bc.ForcingFunctionBC()(x, res, dres);
        v2 = res[0];
    }
    int64_t nvec = datavec.size();

    switch (bc.Type()) {
        case 0 :            // Dirichlet condition
        {
            int64_t firsti = 0;
            for(int i = 0; i<nvec; i++) {
                if(datavec[i].phi.Rows() < 1) continue;
                if(datavec[i].fShapeType != TPZMaterialData::EScalarShape) DebugStop();
                int64_t firstj = 0;
                for(int j = 0; j<nvec; j++) {
                    if(datavec[j].phi.Rows() < 1) continue;
                    ek.AddContribution(firsti,firstj,datavec[i].phi,false,datavec[j].phi,true,weight*fBigNumber);
                    firstj += datavec[j].phi.Rows();
                }
                auto nshape = datavec[i].phi.Rows();
                for(int k = 0; k<nshape; k++) ef(firsti+k,0) += (weight*fBigNumber*v2)*datavec[i].phi(k,0);
                firsti += datavec[i].phi.Rows();
            }
        }
            break;
        case 1 :            // Neumann condition
        {
            int64_t firsti = 0;
            for(int i = 0; i<nvec; i++) {
                auto nshape = datavec[i].phi.Rows();
                if(nshape < 1) continue;
                if(datavec[i].fShapeType != TPZMaterialData::EScalarShape) DebugStop();
                for(int k = 0; k<nshape; k++) ef(firsti+k,0) += (weight*v2)*datavec[i].phi(k,0);
                firsti += datavec[i].phi.Rows();
            }
        }
            break;
        case 2 :        // mixed condition
        {
            int64_t firsti = 0;
            STATE v1 = bc.Val1()(0,0);
            for(int i = 0; i<nvec; i++) {
                int nshapei = datavec[i].phi.Rows();
                if(nshapei < 1) continue;
                if(datavec[i].fShapeType != TPZMaterialData::EScalarShape) DebugStop();
                int64_t firstj = 0;
                for(int j = 0; j<nvec; j++) {
                    int nshapej = datavec[j].phi.Rows();
                    if(nshapej < 1) continue;
                    ek.AddContribution(firsti,firstj,datavec[i].phi,false,datavec[j].phi,true,weight*v1);
                    firstj += datavec[j].phi.Rows();
                }
                auto nshape = datavec[i].phi.Rows();
                for(int k = 0; k<nshape; k++) ef(firsti+k,0) += (weight*v2)*datavec[i].phi(k,0);
                firsti += datavec[i].phi.Rows();
            }
        }
            break;
        default:
            DebugStop();
    }
}

void TPZGFemDarcyFlow::Errors(const TPZVec<TPZMaterialDataT<STATE>> &data, TPZVec<REAL> &errors)
{
    if(!fExactSol) return;

    errors.Resize(NEvalErrors());
    errors.Fill(0.0);

    TPZManVector<STATE> u_exact(1);
    TPZFNMatrix<9,STATE> du_exact;


    if(this->fExactSol){

        this->fExactSol(data[1].x,u_exact,du_exact);
    }
    int nvec = data.size();
    STATE pressure = 0.;
    for(int i = 0; i<nvec; i++) {
        if(data[i].fShapeType != TPZMaterialData::EScalarShape) DebugStop();
        pressure += data[i].sol[0][0];
    }

    // errors[0] norm L2 || u ||_l2

    errors[0] = (pressure-u_exact[0])*(pressure-u_exact[0]);//exact error pressure

    // errors[1] Semi norm H1 || grad u ||_l2

    TPZManVector<STATE,3> sol(1),dsol(3,0.);

    TPZFNMatrix<9,REAL> flux(3,0);
    for(int i = 0; i<nvec; i++) {
        if(data[i].fShapeType != TPZMaterialData::EScalarShape) DebugStop();
        TPZFMatrix<REAL> &dsolaxes = data[i].dsol[0];
        TPZFNMatrix<9,REAL> fluxloc(3,0);
        TPZAxesTools<REAL>::Axes2XYZ(dsolaxes, fluxloc, data[i].axes);
        flux += fluxloc;
    }

    for(int id=0; id<fDim; id++) {
        REAL diff = fabs(flux(id,0) - du_exact(id,0));
        errors[1]  += diff*diff;
    }

    // error[2] H1 norm

    errors[2] = errors[0] +errors[1];

    // error[3] Energy norm || u ||_e = a(u,u)= int_K K gradu.gradu dx

    STATE KPerm = GetPermeability(data[0].x);


    TPZFNMatrix<9,REAL> gradpressure(fDim,1),Kgradu(fDim,1);
    for (int i=0; i<fDim; i++) {
        gradpressure(i,0) = du_exact(i,0);
        Kgradu(i,0) = gradpressure(0)*KPerm;
    }

    REAL energy = 0.;
    for (int i=0; i<fDim; i++) {
        for (int j=0; j<fDim; j++) {
            double cperm =0.;
            if(i==j)
                cperm = KPerm;
            energy += cperm*fabs(flux(j,0) - du_exact(j,0))*fabs(flux(i,0) - du_exact(i,0));
        }
    }

    errors[3] = energy;
}
/**@}*/
