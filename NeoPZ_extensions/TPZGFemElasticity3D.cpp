#include "TPZGFemElasticity3D.h"
#include "pzaxestools.h"


void TPZGFemElasticity3D::Contribute(const TPZVec<TPZMaterialDataT<STATE>> &data,
                                 REAL weight,
                                 TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef) {
    TPZMaterialDataT<STATE> datalocal = data[0];
    datalocal.fShapeType = TPZMaterialData::EVecShape;
    datalocal.sol[0].Fill(0.);
    datalocal.dsol[0].Zero();
    int nshape = ek.Rows();
    datalocal.phi.Redim(3,nshape);
    datalocal.dphix.Redim(9,nshape);
    int shapecount = 0;
    int nvec = data.size();
    for(int ivec = 0; ivec < nvec; ivec++)
    {
        if(datalocal.fNeedsSol) {
            if(data[ivec].fShapeType == TPZMaterialData::EScalarShape || data[ivec].fShapeType == TPZMaterialData::EVecShape) {
                for(int i = 0; i<3; i++) {
                    datalocal.sol[0][i] += data[ivec].sol[0][i];
                    for(int j = 0; j<3; j++) {
                        datalocal.dsol[0](i,j) += data[ivec].dsol[0](i,j);
                    }
                }
            }
        }
        if(data[ivec].fShapeType == TPZMaterialData::EScalarShape) {
            int nshapescalar = data[ivec].phi.Rows();
            for(int ishape = 0; ishape < nshapescalar; ishape++)
            {
                for(int i = 0; i < 3; i++)
                {
                    datalocal.phi(i,shapecount+i) = data[ivec].phi(ishape,0);
                    for(int j = 0; j < 3; j++)
                    {
                        datalocal.dphix(i*3+j,shapecount+i) = data[ivec].dphix(j,ishape);
                    }
                }
                shapecount+=3;
            }
        } else if(data[ivec].fShapeType == TPZMaterialData::EVecShape) {
            int nshapevec = data[ivec].phi.Cols();
            for(int ishape = 0; ishape < nshapevec; ishape++)
            {
                for(int i = 0; i < 3; i++)
                {
                    datalocal.phi(i,shapecount) = data[ivec].phi(i,ishape);
                    for(int j = 0; j < 3; j++)
                    {
                        datalocal.dphix(i*3+j,shapecount) = data[ivec].dphix(i*3+j,ishape);
                    }
                }
                shapecount++;
            }
        }
    }
    if(shapecount != nshape)
    {
        DebugStop();
    }
    TPZElasticity3D::Contribute(datalocal,weight,ek,ef);
}

void TPZGFemElasticity3D::ContributeBC(const TPZVec<TPZMaterialDataT<STATE>> &data,
                                   REAL weight,
                                   TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,
                                   TPZBndCondT<STATE> &bc) {
    TPZMaterialDataT<STATE> datalocal = data[0];
    datalocal.fShapeType = TPZMaterialData::EVecShape;
    datalocal.sol[0].Fill(0.);
    datalocal.dsol[0].Zero();
    int nshape = ek.Rows();
    datalocal.phi.Redim(3,nshape);
    datalocal.dphix.Redim(9,nshape);
    int shapecount = 0;
    int nvec = data.size();
    for(int ivec = 0; ivec < nvec; ivec++)
    {
        if(datalocal.fNeedsSol) {
            if(data[ivec].fShapeType == TPZMaterialData::EScalarShape || data[ivec].fShapeType == TPZMaterialData::EVecShape) {
                for(int i = 0; i<3; i++) {
                    datalocal.sol[0][i] += data[ivec].sol[0][i];
                }
            }
        }
        if(data[ivec].fShapeType == TPZMaterialData::EScalarShape) {
            int nshapescalar = data[ivec].phi.Rows();
            for(int ishape = 0; ishape < nshapescalar; ishape++)
            {
                for(int i = 0; i < 3; i++)
                {
                    datalocal.phi(i,shapecount+i) = data[ivec].phi(ishape,0);
                }
                shapecount+=3;
            }
        } else if(data[ivec].fShapeType == TPZMaterialData::EVecShape) {
            int nshapevec = data[ivec].phi.Cols();
            for(int ishape = 0; ishape < nshapevec; ishape++)
            {
                for(int i = 0; i < 3; i++)
                {
                    datalocal.phi(i,shapecount) = data[ivec].phi(i,ishape);
                }
                shapecount++;
            }
        }
    }
    if(shapecount != nshape)
    {
        DebugStop();
    }
    TPZElasticity3D::ContributeBC(datalocal,weight,ek,ef,bc);
}

void TPZGFemElasticity3D::FillDataRequirements(TPZVec<TPZMaterialDataT<STATE>> &data) const
{
    for(int iv = 0; iv < data.size(); iv++) {
        TPZElasticity3D::FillDataRequirements(data[iv]);
    }
}

void TPZGFemElasticity3D::FillBoundaryConditionDataRequirements(int type,
                                                           TPZVec<TPZMaterialDataT<STATE>> &data) const
{
    for(int iv=0; iv<data.size(); iv++) {
        TPZElasticity3D::FillBoundaryConditionDataRequirements(type, data[iv]);
    }
}




void TPZGFemElasticity3D::Solution(const TPZVec<TPZMaterialDataT<STATE>> &data,
                               int var, TPZVec<STATE> &Solout)
{
    Solout.Fill(0.);
    for(int iv = 0; iv < data.size(); iv++) {
        if(data[iv].fShapeType == TPZMaterialData::EScalarShape || data[iv].fShapeType == TPZMaterialData::EVecShape) {
            TPZManVector<STATE,3> sol(Solout.size(),0.);
            TPZElasticity3D::Solution(data[iv], var, sol);
            for(int i = 0; i < sol.size(); i++) {
                Solout[i] += sol[i];
            }
        }
    }
}


TPZMaterial* TPZGFemElasticity3D::NewMaterial() const
{
    return new TPZGFemElasticity3D(*this);
}


int TPZGFemElasticity3D::ClassId() const
{
    return Hash("TPZGFemElasticity3D") ^
        TBase::ClassId() << 1;
}

void TPZGFemElasticity3D::Read(TPZStream &buf, void *context)
{
    TPZElasticity3D::Read(buf,context);
}
    
void TPZGFemElasticity3D::Write(TPZStream &buf, int withclassid) const
{
    TPZElasticity3D::Write(buf,withclassid);
}


void TPZGFemElasticity3D::Errors(const TPZVec<TPZMaterialDataT<STATE>> &data,
                             TPZVec<REAL> &values) {
    int nvec = data.size();
    int firstvec = 0;
    for(int ivec = 0; ivec < nvec; ivec++)
    {
        if(data[ivec].fShapeType == TPZMaterialData::EScalarShape || data[ivec].fShapeType == TPZMaterialData::EVecShape) {
            firstvec = ivec;
            break;
        }
    }
    for(int ivec = firstvec+1; ivec < nvec; ivec++)
    {

        if(data[ivec].fShapeType == TPZMaterialData::EScalarShape || data[ivec].fShapeType == TPZMaterialData::EVecShape) {
            for(int i = 0; i<3; i++) {
                data[firstvec].sol[0][i] += data[ivec].sol[0][i];
                for(int j = 0; j<3; j++) {
                    data[firstvec].dsol[0](i,j) += data[ivec].dsol[0](i,j);
                }
            }
        }
    }
    TPZElasticity3D::Errors(data[firstvec],values);
}

