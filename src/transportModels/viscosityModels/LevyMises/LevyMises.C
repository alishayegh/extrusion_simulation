/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
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
    viscosityModel(name, viscosityProperties, U, phi),
    rho_("rho", dimMass/pow(dimLength, 3),
    viscosityProperties_),
    sigmaY_("yieldStress", dimPressure, viscosityProperties_),
    strainRateLimit_("strainRateLowerLimit",
        dimless/dimTime, viscosityProperties_),
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
        1./3. * sigmaY_ / (rho_ * strainRateLimit_)
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
//
//    nu_.primitiveFieldRef() = strainRateFunction_->value(sigma());
//
//    volScalarField::Boundary& nuBf = nu_.boundaryFieldRef();
//    const volScalarField::Boundary& sigmaBf = sigma.boundaryField();
//
//    forAll(nuBf, patchi)
//    {
//        nuBf[patchi] = strainRateFunction_->value(sigmaBf[patchi]);
//    }
}

// ************************************************************************* //
