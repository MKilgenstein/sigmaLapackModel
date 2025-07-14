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

#include "LES_CelikEtaIndex.H"
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
    defineTypeNameAndDebug(CelikEtaIndex, 0);
    addToRunTimeSelectionTable(functionObject, CelikEtaIndex, dictionary);
}
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Fcleaunctions  * * * * * * * * * * * //
Foam::tmp<Foam::volScalarField>
Foam::functionObjects::CelikEtaIndex::Delta() const
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
Foam::functionObjects::CelikEtaIndex::eta() const
{
    // Turbulence model, local grid and filter size:
    const momentumTransportModel& turbModel =
            mesh_.lookupType<momentumTransportModel>();

    tmp<volScalarField> tnu(turbModel.nu());
    const volScalarField& nu = tnu();

    const dimensionedScalar epsilonMin(dimArea/pow(dimTime, 3.0), small);

    // (CKJ:Eq. 23)
    return pow025(pow3(nu)/max(epsilon(nu), epsilonMin));
}


Foam::tmp<Foam::volScalarField>
Foam::functionObjects::CelikEtaIndex::epsilon
(
    const volScalarField& nu
) const
{
    const momentumTransportModel& turbModel =
            mesh_.lookupType<momentumTransportModel>();

    tmp<volScalarField> tDelta = Delta();
    const volScalarField& Delta = tDelta();

    tmp<volScalarField> tkSgs(turbModel.k());

    // (Derived based on CKJ:Eq. 25-26, p.031102-5)
    return nuEff(nu, Delta)*tkSgs/(Ck_*sqr(Delta));
}


Foam::tmp<Foam::volScalarField>
Foam::functionObjects::CelikEtaIndex::nuEff
(
    const volScalarField& nu,
    const volScalarField& Delta
) const
{
    const momentumTransportModel& turbModel =
            mesh_.lookupType<momentumTransportModel>();

    tmp<volScalarField> tnut(turbModel.nut());

    // (CKJ:p. 031102-3)
    return nuNum(Delta) + tnut() + nu;
}


Foam::tmp<Foam::volScalarField>
Foam::functionObjects::CelikEtaIndex::nuNum
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
Foam::functionObjects::CelikEtaIndex::kNum
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

Foam::functionObjects::CelikEtaIndex::CelikEtaIndex
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    runTime_(runTime),
    alphaEta_(),
    m_(),
    Cnu_(),
    Cn_(),
    Ck_(),
    CelikEtaIndex_
    (
        IOobject
        (
            "CelikEtaIndex",
            runTime_.name(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE,
            false
        ),
        mesh_,
        dimensionedScalar("CelikEtaIndex", dimless, 0.0),
        calculatedFvPatchField<scalar>::typeName
    )
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::CelikEtaIndex::~CelikEtaIndex()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::CelikEtaIndex::read(const dictionary& dict)
{
    Info << "Reading dictionary entries for Celik eta-based Resolution Index-Model:" << endl;
    alphaEta_ = dict.lookupOrDefault<scalar>("alphaEta", 0.05);
    m_ = dict.lookupOrDefault<scalar>("m", 0.5);
    Cnu_ = dict.lookupOrDefault<scalar>("Cnu", 0.1);
    Cn_ = dict.lookupOrDefault<scalar>("Cn", 1.0);
    Ck_ = dict.lookupOrDefault<scalar>("Ck", 0.0376);
    Info << "Using model parameter: alphaEta = " <<  alphaEta_
         << " | m = " << m_
         << " | Cnu = " << Cnu_
         << " | Cn = " << Cn_
         << " | Ck = " << Ck_
         << endl;

    return true;
}

bool Foam::functionObjects::CelikEtaIndex::execute()
{
    const word& momentumTransportModelName = momentumTransportModel::typeName;
    if 
    (
        mesh_.foundObject<compressible::LESModel>(momentumTransportModelName) ||
        mesh_.foundObject<incompressible::LESModel>(momentumTransportModelName)
    )
    {
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//      Compute the effective Kolmogorov length scale:
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
        Info << "Calculate kolmogorov length scale based quality index:" << endl;
        tmp<volScalarField> teta = eta();
        const volScalarField& eta = teta();

        const volScalarField::Internal h(cbrt(mesh_.V()));

        CelikEtaIndex_.primitiveFieldRef() = 1.0/(1.0 + alphaEta_*pow(h/eta, m_));
        CelikEtaIndex_.correctBoundaryConditions();

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
            << "Unable to find turbulence model in the "
            << "database. Check if the solver is included for foamPostProcess." << exit(FatalError);

        return false;
    }
}


bool Foam::functionObjects::CelikEtaIndex::write()
{
    Info<< tab << "Min: " << tab << min(CelikEtaIndex_).value() 
        << tab << "Max: " << tab << max(CelikEtaIndex_).value() 
        << tab << "Avg: " << tab << average(CelikEtaIndex_).value() << endl;
    
    CelikEtaIndex_.write();

    return true;
}


// ************************************************************************* //
