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
    Foam::functionObjects::resolutionIndex::LES_PopeIndex

Description
    This function object calculates the single-grid quality index of
    Large-Eddy-Simulation based on ratio of resolved to total kinetic 
    enery:
    	
    	M = k_{res}/k_{total}				Eqn.: 10 
    
    with
    
    	k_{res} = (u'.u')/2
    	
    while velocity fluctuation are computed by
    
    	u' = u-\overline{u}
    	
    The total kinetic energy k_{total} is calculated as
    	
    	k_{total} = k_{sgs} + k_{num} + k_{res}
    
    The calculation of the numerical part of kinetic energy k_{num} 
    bases on
    	
    	k_{num} = C_n \left(\frac{h}{\Delta}\right)^2 k_{sgs}
    	
    
    For different flows it can be adapted due changing the model
    constants.
    
    References:
    \verbatim
        Governing equations (tag:P):
            Pope, S. B. (2000).
            Turbulent flows.
            Cambridge, UK: Cambridge Univ. Press
            DOI:10.1017/CBO9780511840531
            
            Governing equations for the denominator kNum term (tag:CKJ):
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
        type            PopeIndex;
        libs            ("libFieldFunctionObjects.so");
        ...
    }
    \endverbatim
    
    The Dictionary entries for adapt the function object to different
    modes are:
    \table
        Property     | Description               | Required | Default value
        type         | type name: PopeIndex   | yes      |
        Cn_	     | Emp. constant (CKJ. 28)	 | no	    |	1.0
    \endtable

SourceFiles
    LES_PopeIndex.C

\*---------------------------------------------------------------------------*/

#ifndef LES_PopeIndex_H
#define LES_PopeIndex_H

#include "fvMeshFunctionObject.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                       Class PopeIndex Declaration
\*---------------------------------------------------------------------------*/

class PopeIndex
:
    public fvMeshFunctionObject
{
    // Private Data

        //- Const reference to the mesh
        //const fvMesh& mesh_;

        //- Database
        const Time& runTime_;

        //- Empirical constant
        scalar Cn_;

        word UName_;
        word UMeanName_;

        //- Quality Index volScalarField
        volScalarField PopeIndex_;


    // Private Member Functions

        //- Return const reference to the mesh
        const fvMesh& mesh() const noexcept
        {
            return mesh_;
        }

        //- Return numerical turbulent kinetic energy field
        tmp<volScalarField> kNum(const volScalarField& kSgs) const;


        tmp<volScalarField> Delta() const;

public:

    //- Runtime type information
    TypeName("PopeIndex");


    // Constructors

        //- Construct from Time and dictionary
        PopeIndex
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        PopeIndex(const PopeIndex&) = delete;


    //- Destructor
    virtual ~PopeIndex();


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
        void operator=(const PopeIndex&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
