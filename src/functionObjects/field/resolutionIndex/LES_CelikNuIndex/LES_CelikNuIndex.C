/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) YEAR OpenFOAM Foundation
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

#include "LES_CelikNuIndex.H"
#include "fvMesh.H"
#include "ListOps.H"
#include "momentumTransportModel.H"
#include "compressibleMomentumTransportModels.H"
#include "incompressibleMomentumTransportModels.H"
#include "LESModel.H"
#include "zeroGradientFvPatchFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(CelikNuIndex, 0);
    addToRunTimeSelectionTable(functionObject, CelikNuIndex, dictionary);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //
Foam::tmp<Foam::volScalarField>
Foam::functionObjects::CelikNuIndex::Delta() const
{
    // LES Model handling:
    const word& momentumTransportModelName = momentumTransportModel::typeName;

    // Check momentum transport model:
    if (mesh_.foundObject<compressible::LESModel>(momentumTransportModelName))
    {
        const compressible::LESModel& lesModel =
            mesh_.lookupObject<compressible::LESModel>(momentumTransportModelName);    
        return lesModel.delta();
    }
    else if (mesh_.foundObject<incompressible::LESModel>(momentumTransportModelName))
    {
        const incompressible::LESModel& lesModel =
            mesh_.lookupObject<incompressible::LESModel>(momentumTransportModelName);
        return lesModel.delta();
    }
}

Foam::tmp<Foam::volScalarField>
Foam::functionObjects::CelikNuIndex::nuNum
(
    const volScalarField& Delta
) const
{
    tmp<volScalarField> tkNum = kNum(Delta);
    const volScalarField& kNum = tkNum.ref();

    // (CKJ:Eq. 35)
    return Foam::sign(kNum)*Cnu_*Delta*Foam::sqrt(kNum);
}


Foam::tmp<Foam::volScalarField>
Foam::functionObjects::CelikNuIndex::kNum
(
    const volScalarField& Delta
) const
{
    const momentumTransportModel& turbModel =
            mesh_.lookupType<momentumTransportModel>();

    tmp<volScalarField> tk(turbModel.k());
    const volScalarField& kSgs = tk();

    volScalarField h
    (
        IOobject
        (
            "h",
            mesh_.time().name(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimLength, small),
        calculatedFvPatchField<scalar>::typeName
    );

    h.internalFieldRef() = cbrt(mesh_.V());

    return Cn_*sqr(h/Delta)*kSgs;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::CelikNuIndex::CelikNuIndex
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    runTime_(runTime),
    alphaNu_(),
    n_(),
    Cnu_(),
    Cn_(),
    CelikNuIndex_
    (
        IOobject
        (
            "CelikNuIndex",
            runTime_.name(),
            mesh_,
            IOobject::NO_READ, //READ_IF_PRESENT,
            IOobject::AUTO_WRITE,
            false
        ),
        mesh_,
        dimensionedScalar("CelikNuIndex", dimless, 0.0),
        calculatedFvPatchField<scalar>::typeName
    )
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::CelikNuIndex::~CelikNuIndex()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::CelikNuIndex::read(const dictionary& dict)
{
    // (Default values from CKJ:p. 031102-3, 031102-5)
    Info << "Reading dictionary entries for Celik nu-based Resolution Index-Model:" << endl;
    alphaNu_ = dict.lookupOrDefault<scalar>("alphaNu", 0.05);
    n_ = dict.lookupOrDefault<scalar>("n", 0.5);
    Cnu_ = dict.lookupOrDefault<scalar>("Cnu", 0.1);
    Cn_ = dict.lookupOrDefault<scalar>("Cn", 1.0);
    Info << "Using model parameter: alphaNu = " <<  alphaNu_
         << " | n = " << n_
         << " | Cnu = " << Cnu_
         << " | Cn = " << Cn_
         << endl;

    return true;
}


bool Foam::functionObjects::CelikNuIndex::execute()
{
    const word& momentumTransportModelName = momentumTransportModel::typeName;
    if 
    (
        mesh_.foundObject<compressible::LESModel>(momentumTransportModelName) ||
        mesh_.foundObject<incompressible::LESModel>(momentumTransportModelName)
    )
    {
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//      Read turbulence model based quantites:
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
        const momentumTransportModel& turbModel =
                mesh_.lookupType<momentumTransportModel>();

        tmp<volScalarField> tnu(turbModel.nu());
        const volScalarField& nu = tnu();

        tmp<volScalarField> tnut(turbModel.nut());
        const volScalarField& nut = tnut();

        tmp<volScalarField> tDelta = Delta();   // filter size

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//      Compute the effective viscosity:
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
        tmp<volScalarField> tnuEff = nuNum(tDelta) + nut + nu;
        const volScalarField& nuEff = tnuEff(); // Effective viscosity

        // Calculate effectiv viscosity-based resolution-index:
        CelikNuIndex_.primitiveFieldRef() = 1.0/(1.0 + alphaNu_*pow(nuEff/nu, n_));
        CelikNuIndex_.correctBoundaryConditions();

        return true;
    }
    else if 
    (
        mesh_.foundType<compressible::RASModel>(momentumTransportModel::typeName) ||
        mesh_.foundType<incompressible::RASModel>(momentumTransportModel::typeName) ||
        mesh_.foundType<compressible::laminarModel>(momentumTransportModel::typeName) ||
        mesh_.foundType<incompressible::laminarModel>(momentumTransportModel::typeName)
    )
    {
        const word& momentumTransportModelName = momentumTransportModel::typeName;
        FatalErrorInFunction()
            << "CelikNuIndex is not available for "
            << momentumTransportModelName
            << " turbulence models."
            << exit(FatalIOError);

        return false;

    }   
    else
    {
        FatalErrorInFunction
            << "Unable to find LES turbulence model in the "
            << "database. Check if the solver is included for foamPostProcess." 
            << exit(FatalError);
        
        return false;
    }
}

bool Foam::functionObjects::CelikNuIndex::write()
{

    Info << "Calculate viscosity based quality index" << endl;

    Info<< tab << "Min: " << tab << min(CelikNuIndex_).value() 
        << tab << "Max: " << tab << max(CelikNuIndex_).value() 
        << tab << "Avg: " << tab << average(CelikNuIndex_).value() << endl;

    CelikNuIndex_.write();

    return true;

}


// ************************************************************************* //
