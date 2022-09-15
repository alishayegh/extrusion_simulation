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

#include "rigidPlastic.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"
#include "fvcGrad.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace viscosityModels
{
    defineTypeNameAndDebug(rigidPlastic, 0);
    addToRunTimeSelectionTable(viscosityModel, rigidPlastic, dictionary);
}
}

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::viscosityModels::rigidPlastic::calcNu() const
/*
{
    volVectorField tU
    (
        //IOobject
        //(
        //    "tU",
        //    U_.time().timeName(),
        //    U_.db(),
        //    IOobject::NO_READ,
        //    IOobject::NO_WRITE
        //),
        //dimensionedVector
        //(
            "tU", dimLength/dimTime, vector::one
        //)
        //U_.mesh()
        //calcNu()
        //dimensionedScalar("tnu", dimLength/pow(dimTime,2),
        //  1.0)
        //1./3. * yieldStress_.value() / 
        //  (strainRateThreshold_.value() *
        //    rho_.value())
    );
    tmp<volScalarField> tsigma = 
        sqrt(2.0)*mag(symm(fvc::grad(
                                        tu
                                    )
                            )
                        );
    //strainRateThreshold_;// strainRate();
    //const volScalarField& sigma = tsigma();

    //forAll(nu_, celli)
    //{
    //    //const double& strainRateThresholdCopy = 
    //    //    strainRateThreshold_;
    //    if ( sigma[celli] > strainRateThreshold_.value() )
    //        nu_[celli] = 1./3. * yieldStress_.value() 
    //            / (sigma[celli] * rho_.value());
    //    else
    //        nu_[celli] = 1./3. * yieldStress_.value() / 
    //            (strainRateThreshold_.value() *
    //                rho_.value());
    //}
    return //max
    //(
    //    VSMALL,
    //    //dimensionedSc1./3. * 
    //    
        yieldStress_ / (max(
                        //strainRate(),
        				strainRateThreshold_,
                        tsigma
                        //dimensionedScalar("strainRateThreshold",
                        //                  dimless/dimTime,
                        //                  strainRateThreshold_)
                        //dimensionedScalar("VSMALL",
                        //                  dimless/dimTime,
                        //                  VSMALL)
                            ) * rho_
                        );
        //yieldStress_ / (strainRateThreshold_ * rho_);
    //);

    //return nu_;
    //return max
    //(
    //    nuMin_,
    //    min
    //    (
    //        nuMax_,
    //        k_*pow
    //        (
    //            max
    //            (
    //                dimensionedScalar("one", dimTime, 1.0)*strainRate(),
    //                dimensionedScalar("VSMALL", dimless, VSMALL)
    //            ),
    //            n_.value() - scalar(1.0)
    //        )
    //    )
    //);
}*/
{}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::viscosityModels::rigidPlastic::rigidPlastic
(
    const word& name,
    const dictionary& viscosityProperties,
    const volVectorField& U,
    const surfaceScalarField& phi
)
:
    viscosityModel(name, viscosityProperties, U, phi),
    rigidPlasticCoeffs_(viscosityProperties.subDict(typeName + "Coeffs")),
    rho_(rigidPlasticCoeffs_.lookup("rho")),
    yieldStress_(rigidPlasticCoeffs_.lookup("yieldStress")),
    strainRateThreshold_(rigidPlasticCoeffs_.lookup("strainRateThreshold")),
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
        //calcNu()
        dimensionedScalar("tnu", dimLength/pow(dimTime,2),
          1.0)
        //1./3. * yieldStress_.value() / 
        //  (strainRateThreshold_.value() *
        //    rho_.value())
    )
{
    correct();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::viscosityModels::rigidPlastic::read
(
    const dictionary& viscosityProperties
)
{
    viscosityModel::read(viscosityProperties);

    viscosityProperties_.lookup("rho") >> rho_;
    viscosityProperties_.lookup("yieldStress") >> yieldStress_;
    viscosityProperties_.lookup("strainRateThreshold") >> 
      strainRateThreshold_;

    return true;
}

//namespace Foam
//{
//namespace viscosityModels
//{

void Foam::viscosityModels::rigidPlastic::correct()
//void rigidPlastic::correct()
{
    tmp<volScalarField> tsigma = strainRate();
    const volScalarField& sigma = tsigma();

    //forAll(rigidPlastic::nu_, celli)
    //{
    //    //const double& strainRateThresholdCopy = 
    //    //    strainRateThreshold_;
    //    if ( sigma[celli] > rigidPlastic::strainRateThreshold_.value() )
    //        rigidPlastic::nu_[celli] = 1./3. * 
    //          rigidPlastic::yieldStress_.value() 
    //            / (sigma[celli] * 
    //              rigidPlastic::rho_.value());
    //    else
    //        nu_[celli] = 1./3. * 
    //          rigidPlastic::yieldStress_.value() / 
    //            (rigidPlastic::strainRateThreshold_.value() *
    //              rigidPlastic::rho_.value());
    //}
    forAll(nu_, celli)
    {
        //const double& strainRateThresholdCopy = 
        //    strainRateThreshold_;
        if ( sigma[celli] > strainRateThreshold_.value() )
            nu_[celli] = 1./3. * 
              yieldStress_.value() 
                / (sigma[celli] * 
                  rho_.value());
        else
            nu_[celli] = 1./3. * 
              yieldStress_.value() / 
                (strainRateThreshold_.value() *
                  rho_.value());
    }
}

//}
//}
// ************************************************************************* //
