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

Class
    Foam::functionObjects::resolutionIndex::LES_CelikEtaIndex

Description
    Single-grid quality index estimator for Large-Eddy-Simulation
    based the cell size to Kolmogorov-Scales ratio :
    	
    	IQ_eta=1/(1+alpha_eta(h/eta)^m)		Eqn.: 10 
    
    For different flows it can be adapted due changing the model
    constants.
    
    References:
    \verbatim
        Governing equations:
            Celik, I., Klein, M., & Janicka, J. (2009).
            Assessment measures for engineering LES applications.
            Journal of fluids engineering, 131(3).
            DOI:10.1115/1.3059703
    \endverbatim

Usage
    Example of function object specification:
    \verbatim
    CelikNuIndex
    {
        type            CelikEtaIndex;
        libs            ("libFieldFunctionObjects.so");
        ...
        
        // Optional entries
        alphaEta    	<scalar>;
        m           	<scalar>;
        Cnu         	<scalar>;
        Cn          	<scalar>;
        Ck          	<scalar>;

    }
    \endverbatim
    
    where the entries mean:
    \table
        Property     | Description               | Required | Default value
        type         | type name: CelikEtaIndex  | yes      |
        alphaEta     | Emp. constant (Eqn. 10)	 | no       |	0.05
        m            | Emp. exponent (Eqn. 10)	 | no       |	0.53
        Ck 	     | Emp. constant (Eqn. 27)	 | no	    |	0.0376
        Cn 	     | Emp. constant (Eqn. 28)	 | no	    |	1.0
        Cnu          | Emp. constant (Eqn. 35)	 | no	    |	0.1
    \endtable

SourceFiles
    LES_CelikEtaIndex.C

\*---------------------------------------------------------------------------*/

#ifndef LES_CelikEtaIndex_H
#define LES_CelikEtaIndex_H

#include "fvMeshFunctionObject.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                       Class CelikEtaIndex Declaration
\*---------------------------------------------------------------------------*/

class CelikEtaIndex
:
    public fvMeshFunctionObject
{
    // Private Data

        //- Const reference to the mesh
        //const fvMesh& mesh_;

        //- Database
        const Time& runTime_;

        //- Empirical constant
        scalar alphaEta_;

        //- Empirical exponent
        scalar m_;

        //- Empirical constant
        scalar Cnu_;

        //- Empirical constant
        scalar Cn_;

        //- Empirical constant
        scalar Ck_;

        //- Quality Index volScalarField
        volScalarField CelikEtaIndex_;


    // Private Member Functions

        //- Return const reference to the mesh
        const fvMesh& mesh() const noexcept
        {
            return mesh_;
        }

        //- Return Kolmogorov length scale field
        tmp<volScalarField> eta() const;

        //- Return kinetic energy dissipation rate field
        tmp<volScalarField> epsilon
        (
            const volScalarField& nu
        ) const;

        //- Return effective eddy viscosity field
        tmp<volScalarField> nuEff
        (
            const volScalarField& nu,
            const volScalarField& Delta
        ) const;

        //- Return numerical eddy viscosity field
        tmp<volScalarField> nuNum(const volScalarField& Delta) const;

        //- Return numerical turbulent kinetic energy field
        tmp<volScalarField> kNum(const volScalarField& Delta) const;

        //- Return Delta field:
        tmp<volScalarField> Delta() const;

public:

    //- Runtime type information
    TypeName("CelikEtaIndex");

    // Constructors

        //- Construct from Time and dictionary
        CelikEtaIndex
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        CelikEtaIndex(const CelikEtaIndex&) = delete;


    //- Destructor
    virtual ~CelikEtaIndex();


    // Member Functions

        //- Read the FUNCTIONOBJECT data
        virtual bool read(const dictionary&);

        //- Return the list of fields required
        virtual wordList fields() const
        {
            return wordList::null();
        }

        //- Execute, currently does nothing
        virtual bool execute();

        //- Write the FUNCTIONOBJECT
        virtual bool write();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const CelikEtaIndex&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
