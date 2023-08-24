/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2019-2021 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "SpalartAllmarasRubberBand.H"
#include "fvOptions.H"
#include "bound.H"
#include "wallDist.H"
#include <iostream>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class BasicTurbulenceModel>
tmp<volScalarField> SpalartAllmarasRubberBand<BasicTurbulenceModel>::chi() const
{
    return nuTilda_/this->nu();
}


template<class BasicTurbulenceModel>
tmp<volScalarField> SpalartAllmarasRubberBand<BasicTurbulenceModel>::fv1
(
    const volScalarField& chi
) const
{
    return pow( scalar(1.0) - exp(scalar(-1.0)*chi/kappa_.value()/scalar(17.0)),2);
}

template<class BasicTurbulenceModel>
tmp<volScalarField::Internal> SpalartAllmarasRubberBand<BasicTurbulenceModel>::fv2
(
    const volScalarField::Internal& chi,
    const volScalarField::Internal& fv1
) const
{
    return scalar(1) - chi/(scalar(1) + chi*fv1);
}

template<class BasicTurbulenceModel>
tmp<volScalarField> SpalartAllmarasRubberBand<BasicTurbulenceModel>::ft2
(
    const volScalarField& chi
) const
{
    const volScalarField chi2(pow(chi,2));

    return Ct3_*exp(-1.0*Ct4_*chi2);
}


template<class BasicTurbulenceModel>
tmp<volScalarField::Internal> SpalartAllmarasRubberBand<BasicTurbulenceModel>::Stilda()
const
{
    const volScalarField chi(this->chi());

    const volScalarField fv1(this->fv1(chi));

    const volScalarField::Internal Omega
    (
        ::sqrt(scalar(2))*mag(skew(fvc::grad(this->U_)().v()))
    );

    return
    (
        Omega + fv2(chi(), fv1())*nuTilda_()/sqr(kappa_*y_)
    );
}


template<class BasicTurbulenceModel>
tmp<volScalarField::Internal> SpalartAllmarasRubberBand<BasicTurbulenceModel>::fw
(
    const volScalarField::Internal& Stilda
) 
{
    Info << "begin fw"<< nl;

    const volScalarField::Internal r
    (
     	max
		(
	 	min
        (
            nuTilda_
           /(
               max
               (
                   Stilda,
                   dimensionedScalar(Stilda.dimensions(), 1e-10)
               )
              *sqr(kappa_*y_)
            ),
            scalar(1000)
        ),
		1e-10
		)
    );

    const volScalarField::Internal g(r + Cw2_*(pow6(r) - r));

    const volScalarField::Internal r2(pow(r,2));

    const volScalarField::Internal F0(
        0.540*r - 0.130*r2
    );

    const volScalarField::Internal F1(
        -0.213 + 0.623*r
    );

    const volScalarField::Internal F2(
        0.049 - 0.365*r + 0.316*r2
    );

    double C1;
    double C2;
    double C3;
    double dC1;
    double dC2;

    double cutX;
    double F0X;
    double F1X;
    double F2X;
    double FwX;
    double dF0;
    double dF1;
    double dF2;
    double dFw;
    

    //Info << "begin fw1.0"<< nl;
    
	C1 = Cb1_.value()/kappa_.value()/kappa_.value()/Cw1_.value();
    C2 = scalar(1)/sigmaNut_.value()/Cw1_.value();
    C3 = (scalar(1) + Cb2_.value())/sigmaNut_.value()/Cw1_.value();

    //Info << "begin fw1.0.1"<< nl;

	const volScalarField::Internal FwUp(C1/max(r,1e-10)+C2*F2/F0/F0+C3*pow(F1/F0,2));

    //Info << "begin fw1.1"<< nl;
   
    cutX = Cw0_.value();

    F0X = 0.540*cutX - 0.130*pow(cutX,2);

    F1X = -0.213 + 0.623*cutX;

    F2X = 0.049 - 0.365*cutX + 0.316*pow(cutX,2);

    //Info << "begin fw1.2"<< nl;
    
	FwX = C1/cutX + C2*F2X/F0X/F0X + C3*pow(F1X/F0X,2);

    //Info << "begin fw1.3"<< nl;

	//C0 continue
	dC1 = FwX/cutX;
	dC2 = 0.0;

    Info << "dC1 = "<<dC1<<", dC2 = "<<dC2<<nl;
    Info << "cutx = "<<cutX<<", FwX = "<<FwX<<nl;

    const volScalarField::Internal FwDown(dC1*r + dC2*r*r);


    volScalarField::Internal fwv(g);

    forAll(fwv,k)
    {
        
        if (r[k] <= 1.0 && r[k]>=cutX){
            fwv[k] = FwUp[k];
        }
        else if (r[k]>=0.0 && r[k]<cutX)
        {
            fwv[k] = FwDown[k];
        }
        else if (r[k]>1.0)
        {
            fwv[k] = (pow(10,2.0*Cs2_.value()-1.0)-1.0)*tanh((r[k]-1.0)*5.0/pow(10,4.0*Cs1_.value()-1.0)) + 1.0;
        }
        else if (r[k]<0.0)
        {
        }
       this->r_[k] = r[k];
	   this->fw_[k] = fwv[k];
    }

    tmp<volScalarField::Internal> fwvv
    (
        new volScalarField::Internal
        (   
            IOobject
            (   
                "fwvv",
                this->runTime_.constant(),
                this->mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),  
            fwv
        )
    );

    return fwvv;
}

template<class BasicTurbulenceModel>
void SpalartAllmarasRubberBand<BasicTurbulenceModel>::correctNut()
{
    this->nut_ = nuTilda_*this->fv1(this->chi());
    this->nut_.correctBoundaryConditions();
    fv::options::New(this->mesh_).correct(this->nut_);

    BasicTurbulenceModel::correctNut();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
SpalartAllmarasRubberBand<BasicTurbulenceModel>::SpalartAllmarasRubberBand
(
    const alphaField& alpha,
    const rhoField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const transportModel& transport,
    const word& propertiesName,
    const word& type
)
:
    eddyViscosity<RASModel<BasicTurbulenceModel>>
    (
        type,
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        transport,
        propertiesName
    ),

    sigmaNut_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "sigmaNut",
            this->coeffDict_,
            scalar(2)/scalar(3)
        )
    ),
    kappa_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "kappa",
            this->coeffDict_,
            0.41
        )
    ),

    Cb1_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "Cb1",
            this->coeffDict_,
            0.1355
        )
    ),
    Cb2_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "Cb2",
            this->coeffDict_,
            0.622
        )
    ),
    Cw1_(Cb1_/sqr(kappa_) + (scalar(1) + Cb2_)/sigmaNut_),
    Cw2_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "Cw2",
            this->coeffDict_,
            0.3
        )
    ),
    Cw3_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "Cw3",
            this->coeffDict_,
            2.0
        )
    ),
    Cv1_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "Cv1",
            this->coeffDict_,
            7.1
        )
    ),
    Ct3_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "Ct3",
            this->coeffDict_,
            1.2
        )
    ),
    Ct4_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "Ct4",
            this->coeffDict_,
            0.5
        )
    ),
    Cs1_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "Cs1",
            this->coeffDict_,
            0.25
        )
    ),
    Cs2_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "Cs2",
            this->coeffDict_,
            1.225
        )
    ),
    Cw0_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "Cw0",
            this->coeffDict_,
            0.4
        )
    ),
    nuTilda_
    (
        IOobject
        (
            "nuTilda",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),
    r_
    (
        IOobject
        (
            "r",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
		this->mesh_
	),
    fw_
    (
        IOobject
        (
            "fw",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
		this->mesh_
	),
    y_(wallDist::New(this->mesh_).y())
{
    if (type == typeName)
    {
        this->printCoeffs(type);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
bool SpalartAllmarasRubberBand<BasicTurbulenceModel>::read()
{
    if (eddyViscosity<RASModel<BasicTurbulenceModel>>::read())
    {
        sigmaNut_.readIfPresent(this->coeffDict());
        kappa_.readIfPresent(this->coeffDict());

        Cb1_.readIfPresent(this->coeffDict());
        Cb2_.readIfPresent(this->coeffDict());
        Cw1_ = Cb1_/sqr(kappa_) + (scalar(1) + Cb2_)/sigmaNut_;
        Cw2_.readIfPresent(this->coeffDict());
        Cw3_.readIfPresent(this->coeffDict());
        Cv1_.readIfPresent(this->coeffDict());
        Ct3_.readIfPresent(this->coeffDict());
        Ct4_.readIfPresent(this->coeffDict());
        Cs1_.readIfPresent(this->coeffDict());
        Cs2_.readIfPresent(this->coeffDict());
        Cw0_.readIfPresent(this->coeffDict());

        return true;
    }

    return false;
}


template<class BasicTurbulenceModel>
tmp<volScalarField> SpalartAllmarasRubberBand<BasicTurbulenceModel>::DnuTildaEff() const
{
    return tmp<volScalarField>::New
    (
        "DnuTildaEff",
        (nuTilda_ + this->nu())/sigmaNut_
    );
}


template<class BasicTurbulenceModel>
tmp<volScalarField> SpalartAllmarasRubberBand<BasicTurbulenceModel>::k() const
{
    // (B:Eq. 4.50)
    const scalar Cmu = 0.09;

    return tmp<volScalarField>::New
    (
        IOobject
        (
            IOobject::groupName("k", this->alphaRhoPhi_.group()),
            this->runTime_.timeName(),
            this->mesh_
        ),
        cbrt(this->fv1(this->chi()))
        *nuTilda_
        *::sqrt(scalar(2)/Cmu)
        *mag(symm(fvc::grad(this->U_))),
        this->nut_.boundaryField().types()
    );
}


template<class BasicTurbulenceModel>
tmp<volScalarField> SpalartAllmarasRubberBand<BasicTurbulenceModel>::epsilon() const
{
    // (B:Eq. 4.50)
    const scalar Cmu = 0.09;
    const dimensionedScalar nutSMALL(sqr(dimLength)/dimTime, SMALL);

    return tmp<volScalarField>::New
    (
        IOobject
        (
            IOobject::groupName("epsilon", this->alphaRhoPhi_.group()),
            this->runTime_.timeName(),
            this->mesh_
        ),
        pow(this->fv1(this->chi()), 0.5)
        *pow(::sqrt(Cmu)*this->k(), 2)
        /(nuTilda_ + this->nut_ + nutSMALL),
        this->nut_.boundaryField().types()
    );
}


template<class BasicTurbulenceModel>
tmp<volScalarField> SpalartAllmarasRubberBand<BasicTurbulenceModel>::omega() const
{
    // (P:p. 384)
    const scalar betaStar = 0.09;
    const dimensionedScalar k0(sqr(dimLength/dimTime), SMALL);

    return tmp<volScalarField>::New
    (
        IOobject
        (
            IOobject::groupName("omega", this->alphaRhoPhi_.group()),
            this->runTime_.timeName(),
            this->mesh_
        ),
        this->epsilon()/(betaStar*(this->k() + k0)),
        this->nut_.boundaryField().types()
    );
}


template<class BasicTurbulenceModel>
void SpalartAllmarasRubberBand<BasicTurbulenceModel>::correct()
{
    if (!this->turbulence_)
    {
        return;
    }

    {
        // Construct local convenience references
        const alphaField& alpha = this->alpha_;
        const rhoField& rho = this->rho_;
        const surfaceScalarField& alphaRhoPhi = this->alphaRhoPhi_;
        fv::options& fvOptions(fv::options::New(this->mesh_));

        eddyViscosity<RASModel<BasicTurbulenceModel>>::correct();

        const volScalarField::Internal Stilda(this->Stilda());
        const volScalarField::Internal ft2(this->ft2(this->chi()));

        tmp<fvScalarMatrix> nuTildaEqn
        (
            fvm::ddt(alpha, rho, nuTilda_)
          + fvm::div(alphaRhoPhi, nuTilda_)
          - fvm::laplacian(alpha*rho*DnuTildaEff(), nuTilda_)
          - Cb2_/sigmaNut_*alpha*rho*magSqr(fvc::grad(nuTilda_))
         ==
            Cb1_*alpha()*rho()*Stilda*nuTilda_()
          - fvm::Sp( alpha()*rho()*( Cw1_*fw(Stilda)*nuTilda_()/sqr(y_) ), nuTilda_)
          + fvOptions(alpha, rho, nuTilda_)
        );

        nuTildaEqn.ref().relax();
        fvOptions.constrain(nuTildaEqn.ref());
        solve(nuTildaEqn);
        fvOptions.correct(nuTilda_);
        bound(nuTilda_, dimensionedScalar(nuTilda_.dimensions(), Zero));
        nuTilda_.correctBoundaryConditions();
    }

    // Update nut with latest available nuTilda
    correctNut();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam
/*
#include "addToRunTimeSelectionTable.H"
#include "makeTurbulenceModel.H"
#include "RASModel.H"
#include "transportModel.H"
#include "incompressibleTurbulenceModel.H"
#include "IncompressibleTurbulenceModel.H"

namespace Foam
{
		    typedef IncompressibleTurbulenceModel<transportModel> transportModelIncompressibleTurbulenceModel;
			typedef RASModel<transportModelIncompressibleTurbulenceModel> RAStransportModelIncompressibleTurbulenceModel;
}

makeTemplatedTurbulenceModel(transportModelIncompressibleTurbulenceModel, RAS, SpalartAllmarasRubberBand)
*/

// ************************************************************************* //
