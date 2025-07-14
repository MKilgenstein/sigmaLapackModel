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

Class
    Foam::compressible::LESModels::Sigma

Description
    Sigma model, 
    according to Franck Nicoud et.al, 2011, Physics of Fluids, vol. 23, "Using singular values to build a subgrid-scale model for large eddy simulations"
    and Baya Toda Hubert et.al, "A dynamic procedure for advanced subgrid-scale models and wall-bounded flows"

    ATTENTION: D = grad(U)

    Algebraic eddy viscosity SGS model founded on the assumption that
    local equilibrium prevails.
    Thus,
    \verbatim
        B = 2/3*k*I - 2*nuSgs*dev(symm(D))
        Beff = 2/3*k*I - 2*nuEff*dev(symm(D))

    where

        k = cI*delta^2*||dev(symm(D)||^2
	    epsilon = nuEff*||dev(symm(D)||^2 (Diss. Brussies, 2014)
        nuSgs = cDSigma*delta^2*DSigma
        nuEff = nuSgs + nu

    In the dynamic version of the choric  Smagorinsky model
    the coefficients cI and cD are calculated during the simulation,

        cI=<K*m>/<m*m>

    and

        cD=1/2*<L.M>/<M.M>,

    where

        K = 0.5*(F(U.U) - F(U).F(U))
        m = delta^2*(4*||F(dev(symm(D)))||^2 - F(||dev(symm(D))||^2))
        L = dev(F(U*U) - F(U)*F(U))
        M = delta^2*(F(DSigma*dev(symm(D))) - 4*F(DSigma)*F(dev(symm(D))))

    The averaging <...> is over the whole domain, i.e. homogeneous turbulence
    is assumed, containing SVS weighting (cf. paper of Baya Toda Hubert et.al)
    \endverbatim

SourceFiles
    sigmaLapackModel.C

\*---------------------------------------------------------------------------*/

#ifndef sigmaLapackModel_H
#define sigmaLapackModel_H

//#include "SVD.H"
#include "LESModel.H"
#include "LESeddyViscosity.H"
#include "simpleFilter.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace LESModels
{

/*---------------------------------------------------------------------------*\
                           Class sigmaLapackModel Declaration
\*---------------------------------------------------------------------------*/

template<class BasicTurbulenceModel>
class sigmaLapackModel
:
    public LESeddyViscosity<BasicTurbulenceModel>
{
    // Private Member Functions
        //- No copy construct
        sigmaLapackModel(const sigmaLapackModel&) = delete;

        //- No copy assignment
        void operator=(const sigmaLapackModel&) = delete;

protected:
    // Protected data

    // Fields:
	    dimensionedScalar cSigma_;
        volScalarField DSigma_;
        volScalarField k_;
	    volScalarField epsilon_;

    // Filters:
        autoPtr<LESfilter> filterPtr_;
        LESfilter& filter_;
    	
    // Protected Member Functions

        //- Update the SGS eddy-viscosity
        void correctNut
        (
            const volTensorField D,
            const volVectorField& U
        );

        virtual void correctNut();

        //- Calculate coefficients cD, cI from filtering velocity field
        tmp<volScalarField> DSigma(const volTensorField& D) const;
//    	scalar cS(const volTensorField& D, volScalarField& DSigma) const;
    	dimensionedScalar cI(const volTensorField& D) const;


public:

    typedef typename BasicTurbulenceModel::alphaField alphaField;
    typedef typename BasicTurbulenceModel::rhoField rhoField;
    typedef typename BasicTurbulenceModel::transportModel transportModel;

    //- Runtime type information
    TypeName("sigmaLapackModel");

    // Constructors

        //- Construct from components
    /*
        sigmaLapackModel
        (
            const volScalarField& rho,
            const volVectorField& U,
            const surfaceScalarField& phi,
            const basicThermo& thermoPhysicalModel
            //transportModel& transport,
            //const word& turbulenceModelName = turbulenceModel::typeName,
            //const word& modelName = typeName
        );
        */
        sigmaLapackModel
        (
            const alphaField& alpha,
            const rhoField& rho,
            const volVectorField& U,
            const surfaceScalarField& alphaRhoPhi,
            const surfaceScalarField& phi,
            const transportModel& transport,
            const word& propertiesName = turbulenceModel::propertiesName,
            const word& type = typeName

        );


    //- Destructor
    virtual ~sigmaLapackModel() = default;


    // Member Functions

        //- Read LESProperties dictionary
        virtual bool read();

        //- Return SGS kinetic energy
        virtual tmp<volScalarField> k() const
        {
            return k_;
        }
        
        //- Return turbulence dissipation rate 
        virtual tmp<volScalarField> epsilon() const
        {
            return epsilon_;
        }

        //- Return the effective diffusivity for k
        tmp<volScalarField> DkEff() const
        {
            return tmp<volScalarField>
            (
                new volScalarField("DkEff", this->nut_ + this->nu())
            );
        }

        //- Correct Eddy-Viscosity and related properties
        virtual void correct();

        
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "sigmaLapackModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
