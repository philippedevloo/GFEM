#include "TPZGFemElasticity2D.h"
#include "pzaxestools.h"
#include <fstream>

void TPZGFemElasticity2D::Contribute(const TPZVec<TPZMaterialDataT<STATE>> &data,
                                 REAL weight,
                                 TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef) {
    TPZFNMatrix<9,STATE> DMatrix(3,3,0.);
    ComputeDMatrix(data[0], DMatrix);
    if(0)
    {
        TPZElasticity2D::Contribute(data[0],weight,ek,ef);
        std::ofstream out("ContributePrev.nb");
        ek.Print("ek = ",out,EMathematicaInput);
        DMatrix.Print("DMatrix = ",out,EMathematicaInput);
        exit(-1);
        return;
    }
    int ndata = data.size();
    int nfunc = 0;
    for(int idata = 0; idata < ndata; idata++) {
        if(data[idata].fShapeType == TPZMaterialData::EScalarShape) {
            nfunc += 2*data[idata].phi.Rows();    
        }
        else if(data[idata].fShapeType == TPZMaterialData::EVecShape) {
            nfunc += data[idata].phi.Cols();
        }
    }
    if(nfunc != ek.Rows()) DebugStop();
    TPZManVector<STATE,3> force(ff);	
	if(fForcingFunction) {            // phi(in, 0) :  node in associated forcing function
		TPZManVector<STATE> res(3);
		fForcingFunction(data[0].x,res);
		force[0] = res[0];
		force[1] = res[1];
		force[2] = res[2];
	}

    TPZFMatrix<REAL> B(3,nfunc,0),DB(3,nfunc,0);
    nfunc = 0;
    bool trigger = false;
    for(int idata = 0; idata < ndata; idata++) {
        if(data[idata].fShapeType == TPZMaterialData::EScalarShape) {
            int nshape = data[idata].phi.Rows();
            TPZFMatrix<REAL> dphixyz(3,nshape);
            TPZAxesTools<REAL>::Axes2XYZ(data[idata].dphix, dphixyz, data[idata].axes);
            for(int i=0; i<nshape; i++) {
                B(0,2*i+nfunc) = dphixyz(0,i);
                B(1,2*i+1+nfunc) = dphixyz(1,i);
                B(2,2*i+nfunc) = dphixyz(1,i);
                B(2,2*i+1+nfunc) = dphixyz(0,i);
                ef(2*i+nfunc,0) += weight*data[idata].phi(i,0)*force[0];
                ef(2*i+1+nfunc,0) += weight*data[idata].phi(i,0)*force[1];
            }
            nfunc += 2*data[idata].phi.Rows();
        }
        else if(data[idata].fShapeType == TPZMaterialData::EVecShape) {
            int nshape = data[idata].phi.Cols();
            TPZFMatrix<REAL> dphixyz(6,nshape);
            VecAxesToXYZ(2, data[idata].axes, data[idata].dphix, dphixyz);
            for(int i=0; i<nshape; i++) {
                B(0,i+nfunc) = dphixyz(0,i);
                B(1,i+nfunc) = dphixyz(4,i);
                B(2,i+nfunc) = (dphixyz(1,i) + dphixyz(3,i));
                ef(i+nfunc,0) += weight*(data[idata].phi(0,i)*force[0]+data[idata].phi(1,i)*force[1]);
            }
            // trigger = true;
            nfunc += data[idata].phi.Cols();
        }
    }
    DB = DMatrix*B;
    ek.AddContribution(0,0,B,true,DB,false,weight);
    if(trigger)
    {
        std::ofstream out("ContributeNew.nb");
        ek.Print("eknew = ",out,EMathematicaInput);
        DMatrix.Print("DMatrix = ",out,EMathematicaInput);
        B.Print("B = ",out,EMathematicaInput);
        DB.Print("DB = ",out,EMathematicaInput);
        out << "Detjac = " << data[0].detjac << std::endl;
        out << "Weight = " << weight << std::endl;
        exit(-1);
        return;
    }

}

void TPZGFemElasticity2D::ContributeBC(const TPZVec<TPZMaterialDataT<STATE>> &data,
                                   REAL weight,
                                   TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,
                                   TPZBndCondT<STATE> &bc) {
    // return TPZElasticity2D::ContributeBC(data[0],weight,ek,ef,bc);
    int nvec = data.size();
    int nshape = ek.Rows();
    int nstate = 2;
    int firstshape = 0;
    TPZFMatrix<REAL> phikeep = data[0].phi;
    TPZFMatrix<REAL> phiglob(nstate,nshape);
    TPZMaterialData::MShapeFunctionType shapetypekeep = data[0].fShapeType;
    for(int iv = 0; iv < nvec; iv++) {
        if(data[iv].fShapeType == TPZMaterialData::EScalarShape) {
            int nphi = data[iv].phi.Rows();
            for(int i=0; i<nphi; i++) {
                phiglob(0,firstshape+2*i) = data[iv].phi(i,0);
                phiglob(1,firstshape+2*i+1) = data[iv].phi(i,0);
            }
            firstshape += 2*nphi;
        } else if (data[iv].fShapeType == TPZMaterialData::EVecShape) {
            int nphi = data[iv].phi.Cols();
            for(int i=0; i<nphi; i++) {
                phiglob(0,firstshape+i) = data[iv].phi(0,i);
                phiglob(1,firstshape+i) = data[iv].phi(1,i);
            }
            firstshape += nphi;
        }
    }
    if(firstshape != nshape) DebugStop();
    // data[0].phi = phiglob;
    // data[0].fShapeType = TPZMaterialData::EVecShape;
    TPZElasticity2D::ContributeBC(data[0],weight,ek,ef,bc);
    data[0].phi = phikeep;
    data[0].fShapeType = shapetypekeep;
}

void TPZGFemElasticity2D::FillDataRequirements(TPZVec<TPZMaterialDataT<STATE>> &data) const
{
    TPZElasticity2D::FillDataRequirements(data[0]);
}

void TPZGFemElasticity2D::FillBoundaryConditionDataRequirements(int type,
                                                           TPZVec<TPZMaterialDataT<STATE>> &data) const
{
    TPZElasticity2D::FillBoundaryConditionDataRequirements(type, data[0]);
}




void TPZGFemElasticity2D::Solution(const TPZVec<TPZMaterialDataT<STATE>> &data,
                               int var, TPZVec<STATE> &Solout)
{
    int nvec = data.size();
    TPZVec<STATE> Soloutloc(Solout.size(),0.);
    Solout.Fill(0.);
    for(int ivec = 0; ivec < nvec; ivec++) {
        if(data[ivec].fShapeType != TPZMaterialData::EEmpty) {
            TPZElasticity2D::Solution(data[ivec], var, Soloutloc);
            for(int i=0; i<Soloutloc.size(); i++) Solout[i] += Soloutloc[i];
        }
    }
}


TPZMaterial* TPZGFemElasticity2D::NewMaterial() const
{
    return new TPZGFemElasticity2D(*this);
}


int TPZGFemElasticity2D::ClassId() const
{
    return Hash("TPZGFemElasticity2D") ^
        TPZElasticity2D::ClassId() << 1;
}

void TPZGFemElasticity2D::Read(TPZStream &buf, void *context)
{
    TPZElasticity2D::Read(buf,context);
}
	
void TPZGFemElasticity2D::Write(TPZStream &buf, int withclassid) const
{
    TPZElasticity2D::Write(buf,withclassid);
}


void TPZGFemElasticity2D::Errors(const TPZVec<TPZMaterialDataT<STATE>> &data,
                             TPZVec<REAL> &values) {
    TPZElasticity2D::Errors(data[1],values);
}

#include "pzaxestools.h"

    // transform derivatives in the direction of axes to derivatives in xy
void TPZGFemElasticity2D::VecAxesToXYZ(int nstate, TPZFMatrix<REAL> &axes, TPZFMatrix<REAL> &dudaxes, TPZFMatrix<REAL> &dudxyz)
{
    int dim = axes.Rows();
    if(dudaxes.Rows() != nstate*dim) DebugStop();
    if(dudxyz.Rows() != nstate*3) DebugStop();
    if(dudaxes.Cols() != dudxyz.Cols()) DebugStop();
    int cols = dudxyz.Cols();
    TPZFMatrix<REAL> daxes(dim,cols), dudxloc(3,cols);
    for(int d = 0; d<nstate; d++) {
        for(int d2=0; d2<dim; d2++) for(int c=0; c<cols; c++) {
            daxes(d2,c) = dudaxes(d2+d*dim,c);
        }
        TPZAxesTools<REAL>::Axes2XYZ(daxes, dudxloc, axes);
        for(int d2=0; d2<3; d2++) for(int c=0; c<cols; c++) {
            dudxyz(d2+d*3,c)=dudxloc(d2,c);
        }
    }
}

    // transform derivatives in the direction of xy to derivatives in axes
void TPZGFemElasticity2D::VecXYZToAxes(int nstate, TPZFMatrix<REAL> &axes, TPZFMatrix<REAL> &dudxyz, TPZFMatrix<REAL> &dudaxes)
{
    int dim = axes.Rows();
    if(dudaxes.Rows() != nstate*dim) DebugStop();
    if(dudxyz.Rows() != nstate*3) DebugStop();
    if(dudaxes.Cols() != dudxyz.Cols()) DebugStop();
    int cols = dudxyz.Cols();
    TPZFMatrix<REAL> daxes(dim,cols), dudxloc(3,cols);
    for(int d = 0; d<nstate; d++) {
        for(int d2=0; d2<3; d2++) for(int c=0; c<cols; c++) {
            dudxloc(d2,c) = dudxyz(d2+d*3,c);
        }
        TPZAxesTools<REAL>::XYZ2Axes(daxes, dudxloc, axes);
        for(int d2=0; d2<dim; d2++) for(int c=0; c<cols; c++) {
            dudaxes(d2+d*dim,c)=daxes(d2,c);
        }
    }
}

    // transform derivatives in the direction of axes to derivatives in xy
void TPZGFemElasticity2D::VecToAtomic(int nstate, TPZMaterialData &data)
{
    if(data.fShapeType == TPZMaterialData::EScalarShape) return;
    if(data.fShapeType != TPZMaterialData::EVecShape) DebugStop();
    int nshape = data.phi.Cols()/nstate;
    int dim = data.axes.Rows();
    data.phi.Redim(nshape,1);
    data.dphix.Redim(dim,nshape);
    data.fShapeType = TPZMaterialData::EScalarShape;
}

    // transform a scalar shape structure to a vectorial shape structure
void TPZGFemElasticity2D::AtomicToVec(int nstate, TPZMaterialData &data)
{
    if(data.fShapeType == TPZMaterialData::EVecShape) return;
    if(data.fShapeType != TPZMaterialData::EScalarShape) DebugStop();
    int nshape = data.phi.Rows();
    int dim = data.dphix.Rows();
    TPZFMatrix<STATE> dphivec(nstate*dim,nshape*nstate,0.);
    TPZFMatrix<STATE> phivec(nstate,nshape*nstate,0.);
    for(int i=0; i<nshape; i++)
    {
        for(int jst = 0; jst<nstate; jst++)
        {
            phivec(jst,i*nstate+jst) = data.phi(i,0);
            for(int k=0; k<dim; k++)
            {
                dphivec(k+jst*dim,i*nstate+jst) = data.dphix(k,i);
            }
        }
    }
    data.phi = phivec;
    data.dphix = dphivec;
    data.fShapeType = TPZMaterialData::EVecShape;
}

void TPZGFemElasticity2D::ComputeDMatrix(TPZMaterialData &data, TPZFMatrix<STATE> &DMatrix) const
{
        REAL E(fE_def), nu(fnu_def);
    
    if (fElasticity) {
        TPZManVector<STATE,2> result(2);
        TPZFNMatrix<4,STATE> Dres(0,0);
        fElasticity(data.x, result, Dres);
        E = result[0];
        nu = result[1];
    }
    
    REAL Eover1MinNu2 = E/(1-nu*nu);
    REAL Eover21PlusNu = E/(2.*(1+nu));
    REAL EComp = E*(1-nu)/((1-2*nu)*(1+nu));
    DMatrix.Redim(3,3);
    if(this->fPlaneStress != 1) {
        DMatrix(0,0) = EComp;
        DMatrix(0,1) = EComp*nu/(1-nu);
        DMatrix(0,2) = 0.;
        DMatrix(1,0) = EComp*nu/(1-nu);
        DMatrix(1,1) = EComp;
        DMatrix(1,2) = 0.;
        DMatrix(2,0) = 0.;
        DMatrix(2,1) = 0.;
        DMatrix(2,2) = E/(2.*(1+nu));
    } else {
        DMatrix(0,0) = Eover21PlusNu;
        DMatrix(0,1) = Eover21PlusNu*nu;
        DMatrix(0,2) = 0.;
        DMatrix(1,0) = Eover21PlusNu*nu;
        DMatrix(1,1) = Eover21PlusNu;
        DMatrix(1,2) = 0.;
        DMatrix(2,0) = 0.;
        DMatrix(2,1) = 0.;
        DMatrix(2,2) = Eover21PlusNu*(1.-nu)/2.;
    }
}