/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2021 OpenFOAM Foundation
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

#include "makeCompressibleMomentumTransportModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makeBaseMomentumTransportModel
(
    geometricOneField,
    volScalarField,
    compressibleMomentumTransportModel
);


// -------------------------------------------------------------------------- //
// Laminar models
// -------------------------------------------------------------------------- //

/*#include "Stokes.H"
makeLaminarModel(Stokes);

#include "generalisedNewtonian.H"
makeLaminarModel(generalisedNewtonian);

#include "lambdaThixotropic.H"
makeLaminarModel(lambdaThixotropic);

#include "Maxwell.H"
makeLaminarModel(Maxwell);

#include "Giesekus.H"
makeLaminarModel(Giesekus);

#include "PTT.H"
makeLaminarModel(PTT);*/


// -------------------------------------------------------------------------- //
// RAS models
// -------------------------------------------------------------------------- //

/*#include "kEpsilonPhitF.H"
makeRASModel(kEpsilonPhitF);

#include "kEpsilonZetaF.H"
makeRASModel(kEpsilonZetaF);

#include "EBRSM.H"
makeRASModel(EBRSM);*/


// -------------------------------------------------------------------------- //
// LES models
// -------------------------------------------------------------------------- //

/*#include "sigma.H"
makeLESModel(sigma);

#include "dynamicSmagorinsky.H"
makeLESModel(dynamicSmagorinsky);

#include "dynamicSigma.H"
makeLESModel(dynamicSigma);

#include "shearImprovedSmagorinskyModel.H"
makeLESModel(shearImprovedSmagorinskyModel);

#include "CSM.H"
makeLESModel(CSM);*/

#include "sigmaLapackModel.H"
makeLESModel(sigmaLapackModel)

// ************************************************************************* //
