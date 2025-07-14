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

#include "LES_PopeIndex.H"
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
    defineTypeNameAndDebug(PopeIndex, 0);
    addToRunTimeSelectionTable(functionObject, PopeIndex, dictionary);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //
Foam::tmp<Foam::volScalarField>
Foam::functionObjects::PopeIndex::Delta() const
{
    // LES Model handling: 
    typedef compressible::LESModel cmpModel;
    typedef incompressible::LESModel icoModel;

    const word& momentumTransportModelName = momentumTransportModel::typeName;

    // Check if momentum transport model is of type compressible or 
    // incompressible:
    if (mesh_.foundObject<cmpModel>(momentumTransportModelName))
    {
	const cmpModel& lesModel =
	    mesh_.lookupObject<cmpModel>(momentumTransportModelName);    
        return lesModel.delta();
    }
    else if (mesh_.foundObject<icoModel>(momentumTransportModelName))
    {
	const icoModel& lesModel =
            mesh_.lookupObject<icoModel>(momentumTransportModelName);
        return lesModel.delta();
    }
}


Foam::tmp<Foam::volScalarField>
Foam::functionObjects::PopeIndex::kNum
(
    const volScalarField& kSgs
) const
{
    // Read LES-Model and filter sizes:
    tmp<volScalarField> tDelta = Delta();
    const volScalarField& Delta = tDelta();

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

Foam::functionObjects::PopeIndex::PopeIndex
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    runTime_(runTime),
    Cn_(),
    UName_(),
    UMeanName_(),
    PopeIndex_
    (
        IOobject
        (
            "PopeIndex",
            runTime_.name(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE,
            false
        ),
        mesh_,
        dimensionedScalar("PopeIndex", dimless, 0.0),
        calculatedFvPatchField<scalar>::typeName
    )
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::PopeIndex::~PopeIndex()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::PopeIndex::read(const dictionary& dict)
{
    Info << "Reading dictionary entries:" << endl;
    Cn_ = dict.lookupOrDefault<scalar>("Cn", 1.0);
    UName_ = dict.lookupOrDefault<word>("UName", "U");
    UMeanName_ = dict.lookupOrDefault<word>("UMeanName", "UMean");

    return true;
}

bool Foam::functionObjects::PopeIndex::execute()
{
    Info << "Calculate subgrid-activity scale based quality index:" << endl;
    
    // Check if momentum transport model is RANS-based:
    if 
        (
            mesh_.foundType<compressible::RASModel>(momentumTransportModel::typeName) ||
            mesh_.foundType<incompressible::RASModel>(momentumTransportModel::typeName)
        )
        {
            FatalErrorInFunction()
                << "CelikEtaIndex is not available for RANS-based turbulence models."
                << exit(FatalIOError);

            return false;

        }
    
    // Calculate resolved k field
    const volVectorField& U = mesh_.lookupObject<volVectorField>(UName_);
    
    // Read UMean:
    volVectorField UMean
    (
        IOobject
        (
            UMeanName_,
            runTime_.name(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedVector(UMeanName_, dimVelocity, vector::zero)
    );

    // Load momentum transport model:
    const momentumTransportModel& turbModel =
            mesh_.lookupType<momentumTransportModel>();
    
    // Calculate resolved kinetic energy as k_res = (u'.u')/2 with u' = u - \overline{u}
    const volScalarField kRes(0.5*magSqr(U - UMean));

    // Get subgrid-scale k field from momentum transport model:
    tmp<volScalarField> tk(turbModel.k());
    const volScalarField& kSgs = tk();
    
    // Calculate numerical introduced numerical kinetic energy
    tmp < volScalarField> tkNum = kNum(kSgs);
    const volScalarField& kNum = tkNum();

    // Calculate Pope-Index:
    PopeIndex_.primitiveFieldRef() = kRes/max(kRes+kSgs+kNum, dimensionedScalar(sqr(dimLength/dimTime), small));
    PopeIndex_.correctBoundaryConditions();

    return true;
}


bool Foam::functionObjects::PopeIndex::write()
{
    Info<< tab << "Min: " << tab << min(PopeIndex_).value() 
        << tab << "Max: " << tab << max(PopeIndex_).value() 
        << tab << "Avg: " << tab << average(PopeIndex_).value() << endl;
    PopeIndex_.write();

    return true;
}


// ************************************************************************* //
