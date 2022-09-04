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

\*---------------------------------------------------------------------------*/

#include "LevyMises.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace viscosityModels
{
    defineTypeNameAndDebug(LevyMises, 0);
    addToRunTimeSelectionTable(viscosityModel, LevyMises, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::viscosityModels::LevyMises::LevyMises
(
    const word& name,
    const dictionary& viscosityProperties,
    const volVectorField& U,
    const surfaceScalarField& phi
)
:
// Start foam-extend-4.1
//    viscosityModel(name, viscosityProperties, U, phi),
//    nu0_(viscosityProperties_.lookup("nu")),
//    nu_
//    (
//        IOobject
//        (
//            name,
//            U_.time().timeName(),
//            U_.db(),
//            IOobject::NO_READ,
//            IOobject::NO_WRITE
//        ),
//        U_.mesh(),
//        nu0_
//    )
//{}
// End foam-extend-4.1
    viscosityModel(name, viscosityProperties, U, phi),
    rho_(viscosityProperties_.lookup("rho")),
    sigmaY_(viscosityProperties_.lookup("yieldStress")),
    strainRateLimit_(viscosityProperties_.lookup("strainRateLowerLimit")),
    //rho_("rho", dimMass/pow(dimLength, 3),
    //viscosityProperties_),
    //sigmaY_("yieldStress", dimPressure, viscosityProperties_),
    //strainRateLimit_("strainRateLowerLimit",
    //    dimless/dimTime, viscosityProperties_),
    nu_
    (
        IOobject
        (
            name,
            U_.time().timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        U_.mesh(),
        correct()
        //1./3. * sigmaY_ / (rho_ * strainRateLimit_)
        //1
    )
{
    correct();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::viscosityModels::LevyMises::read
(
    const dictionary& viscosityProperties
)
{
    viscosityModel::read(viscosityProperties);
// Start foam-extend-4.1
    //viscosityProperties_.lookup("nu") >> nu0_;
    //nu_ = nu0_;

    //return true;
// End foam-extend-4.1
    viscosityProperties_.lookup("rho") >>
        rho_;
    viscosityProperties_.lookup("yieldStress") >>
        sigmaY_;
    viscosityProperties_.lookup("strainRateLowerLimit")
        >> strainRateLimit_;
    //nu_ = nu0_;

    return true;

    // From strainRateFunction.C
    //viscosityModel::read(viscosityProperties);

    //LevyMisesCoeffs = viscosityProperties.optionalSubDict
    //(
    //    typeName + "Coeffs"
    //);

    //strainRateFunction_.clear();
    //strainRateFunction_ = Function1<scalar>::New
    //(
    //    "function",
    //    strainRateFunctionCoeffs_
    //);

    //return true;
}

void Foam::viscosityModels::LevyMises::correct()
{
    tmp<volScalarField> tsigma = strainRate();
    const volScalarField& sigma = tsigma();

    forAll(nu_, celli)
    {
        //const double& strainRateLimitCopy = 
        //    strainRateLimit_;
        if ( sigma[celli] > strainRateLimit_.value() )
            nu_[celli] = 1./3. * sigmaY_.value() 
                / (sigma[celli] * rho_.value());
        else
            nu_[celli] = 1./3. * sigmaY_.value() / 
                (strainRateLimit_.value() *
                    rho_.value());
    }
    //nu_.correctBoundaryConditions();
}

// ************************************************************************* //
