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
    Foam::functionObjects::resolutionIndex::LES_CelikNuIndex

Description
    This function object calculates the single-grid quality index of
    Large-Eddy-Simulation based the effective viscosity:
    	
    	IQ_\nu=1/(1+alpha_\nu(nu_{Eff}/nu)^n)		            Eqn.: 10 
    
    where the effective viscosity is

        nu_{Eff} = nu + nu_sgs + nu_num                           

    Therefor the numerical introduced viscosity is defined as

        nu_{num} = C_\nu \Delta \sqrt{k_num}                       Eqn: 35

    while the numerical kinteic energy is

        k_{num} = C_n \left(\frac{h}{\Delta}\right)^2 k_{sgs}      Eqn: 28

    For different flows it can be adapted due changing the model
    constants alpha_nu, C_nu and C_n as well as the exponent n.


    References:
    \verbatim
        Governing equations (tag:CKJ):
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
        type            CelikNuIndex;
        libs            ("libFieldFunctionObjects.so");
        ...

        // Optional entries
        alphaNu     <scalar>;
        n           <scalar>;
        Cnu         <scalar>;
        Cn          <scalar>;
    }
    \endverbatim
    
    where the entries mean:
    \table
        Property     | Description               | Required | Default value
        type         | type name: CelikNuIndex   | yes      |
        alphaNu      | Emp. constant (Eqn. 10)	 | no       |	0.05
        n            | Emp. exponent (Eqn. 10)	 | no       |	0.53
        Cn 	     | Emp. constant (Eqn. 28)	 | no	    |	1.0
        Cnu          | Emp. constant (Eqn. 35)	 | no	    |	0.1
    \endtable

SourceFiles
    LES_CelikNuIndex.C

\*---------------------------------------------------------------------------*/

#ifndef CelikNuIndex_H
#define CelikNuIndex_H

#include "fvMeshFunctionObject.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                       Class CelikNuIndex Declaration
\*---------------------------------------------------------------------------*/

class CelikNuIndex
:
    public fvMeshFunctionObject
{
    // Private Data

        //- Database
        const Time& runTime_;

        //- Empirical constant
        scalar alphaNu_;

        //- Empirical exponent
        scalar n_;

        //- Empirical constant
        scalar Cnu_;

        //- Empirical constant
        scalar Cn_;

        //- Quality Index volScalarField
        volScalarField CelikNuIndex_;



    // Private Member Functions

        //- Return const reference to the mesh
        const fvMesh& mesh() const noexcept
        {
            return mesh_;
        }

        //- Return kinetic energy dissipation rate field
        tmp<volScalarField> epsilon() const;

        //- Return effective eddy viscosity field
        tmp<volScalarField> nuEff() const;

        //- Return numerical eddy viscosity field
        tmp<volScalarField> nuNum(const volScalarField& Delta) const;

        //- Return numerical turbulent kinetic energy field
        tmp<volScalarField> kNum(const volScalarField& Delta) const;

        //- Return Delta field:
        tmp<volScalarField> Delta() const;

public:

    //- Runtime type information
    TypeName("CelikNuIndex");


    // Constructors

        //- Construct from Time and dictionary
        CelikNuIndex
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        CelikNuIndex(const CelikNuIndex&) = delete;


    //- Destructor
    virtual ~CelikNuIndex();


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
        void operator=(const CelikNuIndex&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
