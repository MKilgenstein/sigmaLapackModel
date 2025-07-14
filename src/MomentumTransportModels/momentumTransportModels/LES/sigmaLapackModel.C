/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
     \\/     M anipulation  |
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

#include "sigmaLapackModel.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "volMesh.H"
#include "bound.H"
#include </usr/include/lapack.h>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace LESModels
{


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class BasicTurbulenceModel>
tmp<volScalarField> sigmaLapackModel<BasicTurbulenceModel>::DSigma
(
     const volTensorField& D
) const
{
    //Info << "Calculate eigenvalues of the velocity gradient tensor gradU." << endl;
    /*volScalarField Dsigma(mag(gradU));
    scalarField& DsgCells = Dsigma.primitiveFieldRef();*/

    volScalarField Dsigma(mag(fvc::grad(this->U_)));
    scalarField& DsgCells = Dsigma.primitiveFieldRef();

    forAll(D, celli)
    {
        tensor G = D[celli].T()&D[celli];
        
        // LAPACK erfordert ein Array, daher formatieren wir das Tensorfeld G als float-Array
        double G_array[9] = {G.xx(), G.xy(), G.xz(), 
                             G.xy(), G.yy(), G.yz(), 
                             G.xz(), G.yz(), G.zz()};

        Info << "G = [" << G.xx() << ", " << G.xy() << ", " << G.xz() << "; "
                << G.xy() << ", " << G.yy() << ", " << G.yz() << "; "
                << G.xz() << ", " << G.yz() << ", " << G.zz() << "]" << endl;

        int n = 3;                      // Matrixdimension 3x3
        int lda = 3;                    // leading dimension der Matrix G
        double lambda[3];           	// Array für die Eigenwerte
        double work[8];                 // Arbeitsarray, Größe ≥ 3*n - 1
        int lwork = 15;                 // Größe des Arbeitsarrays
        int info;                       // Fehlercode

        // Berechnen Sie die Eigenwerte und Eigenvektoren von A
        dsyev_("N", "U", &n, G_array, &lda, lambda, work, &lwork, &info, 1, 1);

        if (info==0)
        {
            const scalar lambda_1 = lambda[2];
            const scalar lambda_2 = lambda[1];
            const scalar lambda_3 = lambda[0];

            Info << "Calculate sigular values sigma for lambda_1 = " << lambda_1 
                 << " lambda_2 = " << lambda_2 
                 << " lambda_3 = " << lambda_3 << endl;

            scalar sigma_1 = sqrt(max(SMALL, lambda_1));
            scalar sigma_2 = sqrt(max(SMALL, lambda_2));
            scalar sigma_3 = sqrt(max(SMALL, lambda_3));

            DsgCells[celli] = sigma_3*(sigma_1-sigma_2)*(sigma_2-sigma_3)/(sqr(sigma_1)+VSMALL);
            Info << "DSimga[" << celli << "] = " << DsgCells[celli] << endl;
        }
        else
        {
            WarningInFunction
                << "LAPACK ssyev failed for cell " << celli << " with info = " << info << nl
                << "Setting Dsigma to zero for this cell." << endl;
            DsgCells[celli] = 0.0;
        }
    }

    return tmp<volScalarField>(new volScalarField(Dsigma));
}

template<class BasicTurbulenceModel>
dimensionedScalar sigmaLapackModel<BasicTurbulenceModel>::cI
(
   const volTensorField& D
) const
{
    const volScalarField mm
    (
        sqr(this->delta())*(4*sqr(mag(filter_(dev(symm(D))))) - filter_(sqr(mag(dev(symm(D))))))
    );

    dimensionedScalar mmmm = average(magSqr(mm));

    if (mmmm.value() > VSMALL)
    {
        tmp<volScalarField> KK =
            0.5*(filter_(magSqr(this->U())) - magSqr(filter_(this->U())));

        return average(KK*mm)/mmmm;
    }
    else
    {
        return 0.0;
    }
}

template<class BasicTurbulenceModel>
void sigmaLapackModel<BasicTurbulenceModel>::correctNut()
{
    correctNut(fvc::grad(this->U_), this->U_);
}


template<class BasicTurbulenceModel>
void sigmaLapackModel<BasicTurbulenceModel>::correctNut
(
    const volTensorField D,
    const volVectorField& U
)
{   
    DSigma_ = DSigma(D);
    const dimensionedScalar lim(sqr(dimLength)/dimTime, 1e-10);
    this->nut_ = max(sqr(cSigma_*this->delta()) * DSigma_, lim);	//eq. 7 in TSFP_BayaToda
    this->nut_.correctBoundaryConditions();

    BasicTurbulenceModel::correctNut();
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
sigmaLapackModel<BasicTurbulenceModel>::sigmaLapackModel
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
    LESeddyViscosity<BasicTurbulenceModel>
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
    
    DSigma_
    (
    	IOobject
    	(
        	"DSigma_",
		    this->runTime_.timeName(),
            this->mesh_,
        	IOobject::NO_READ,
        	IOobject::AUTO_WRITE
    	),
    	this->mesh_,
	    dimensionedScalar("zero", dimensionSet(0, 0, -1, 0, 0, 0, 0), 0.0)
    ), 

    k_
    (
        IOobject
        (
            "k",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),

    epsilon_
    (
        IOobject
        (
            "epsilon",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_,
	    dimensionedScalar("zero", dimensionSet(0, 2, -3, 0, 0, 0, 0),0.0)
    ),  

    filterPtr_(LESfilter::New(this->mesh_, this->coeffDict_)),
    filter_(filterPtr_())
    
{
    bound(k_,  this->kMin_);

    bound(epsilon_, this->epsilonMin_);

    if (type == typeName)
    {
        this->printCoeffs(type);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
bool sigmaLapackModel<BasicTurbulenceModel>::read()
{
    if (LESeddyViscosity<BasicTurbulenceModel>::read())
    {
        filter_.read(this->coeffDict());
        cSigma_ = this->coeffDict().lookupOrDefault("cSigma", 1.5);
	
        return true;
    }
    else
    {
        return false;
    }
}

template<class BasicTurbulenceModel>
void sigmaLapackModel<BasicTurbulenceModel>::correct()
{
    if (!this->turbulence_)
    {
        return;
    }

    // Local references
    const volVectorField& U = this->U_;
    const volTensorField D(fvc::grad(U));

    LESeddyViscosity<BasicTurbulenceModel>::correct();

    const volScalarField magSqrDevSymmD(magSqr(dev(symm(D))));

    k_ = cI(D)*sqr(this->delta())*magSqrDevSymmD;
    bound(k_,  this->kMin_);

    //Eq (5.65) in Sagaut, LES for incompressible Fluid
    epsilon_ = filter_(2.0*this->nuEff()*magSqrDevSymmD);
    bound(epsilon_, this->epsilonMin_);

    correctNut();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace Foam

// ************************************************************************* //
