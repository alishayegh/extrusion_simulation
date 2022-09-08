/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::viscosityModels::rigidPlastic

Description
    An incompressible rigidPlastic viscosity model.

SourceFiles
    rigidPlastic.C

\*---------------------------------------------------------------------------*/

#ifndef rigidPlastic_H
#define rigidPlastic_H

#include "viscosityModel.H"
#include "dimensionedScalar.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace viscosityModels
{

/*---------------------------------------------------------------------------*\
                           Class rigidPlastic Declaration
\*---------------------------------------------------------------------------*/

class rigidPlastic
:
    public viscosityModel
{
    // Private data

		dictionary rigidPlasticCoeffs_;

        dimensionedScalar rho_;
        dimensionedScalar yieldStress_;
        dimensionedScalar strainRateThreshold_;

        volScalarField nu_;


    // Private Member Functions

        tmp<volScalarField> calcNu() const;


public:

    //- Runtime type information
    TypeName("rigidPlastic");


    // Constructors

        //- construct from components
        rigidPlastic
        (
            const word& name,
            const dictionary& viscosityProperties,
            const volVectorField& U,
            const surfaceScalarField& phi
        );


    // Destructor

        virtual ~rigidPlastic()
        {}


    // Member Functions

        //- Return the laminar viscosity
        virtual const volScalarField& nu() const
        {
            return nu_;
        }

        //- Correct the laminar viscosity (not appropriate, viscosity constant)
        virtual void correct();
        //{
        //    nu_ = calcNu();
        //}

        //- Read transportProperties dictionary
        virtual bool read(const dictionary& viscosityProperties);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace viscosityModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
