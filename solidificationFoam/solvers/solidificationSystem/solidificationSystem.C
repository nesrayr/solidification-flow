/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2018 OpenFOAM Foundation
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

#include "solidificationSystem.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(solidificationSystem, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solidificationSystem::solidificationSystem
(
    const volVectorField& U
)
:
    IOdictionary
    (
        IOobject
        (
            "transportProperties",
            U.time().constant(),
            U.db(),
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),
    twoPhaseMixture(U.mesh(), *this),

    rho1_
    (
        "rho",
        dimensionSet(1, -3, 0, 0, 0),
        subDict(phase1Name_)
    ),
    rho2_
    (
        "rho",
        dimensionSet(1, -3, 0, 0, 0),
        subDict(phase2Name_)
    ),

    Cp1_
    (
        "Cp",
        dimensionSet(0, 2, -2, -1, 0),
        subDict(phase1Name_)
    ),
    Cp2_
    (
        "Cp",
        dimensionSet(0, 2, -2, -1, 0),
        subDict(phase2Name_)
    ),

    kappa1_
    (
        "kappa",
        dimensionSet(1, 1, -3, -1, 0),
        subDict(phase1Name_)
    ),
    kappa2_
    (
        "kappa",
        dimensionSet(1, 1, -3, -1, 0),
        subDict(phase2Name_)
    ),

    mu1_
    (
        "mu",
        dimensionSet(1, -1, -1, 0, 0),
        subDict(phase1Name_)
    ),
    mu2_
    (
        "mu",
        dimensionSet(1, -1, -1, 0, 0),
        subDict(phase2Name_)
    ),

    D1_
    (
        "D",
        dimensionSet(0, 2, -1, 0, 0),
        subDict(phase1Name_)
    ),
    D2_
    (
        "D",
        dimensionSet(0, 2, -1, 0, 0),
        subDict(phase2Name_)
    ),
    DAS_
    (
        "DAS",
        dimensionSet(0, 1, 0, 0, 0),
        subDict(phase1Name_)
    ),
    betaT_
    (
        "betaT",
        dimensionSet(0, 0, 0, -1, 0),
        subDict(phase2Name_)
    ),
    betaC_
    (
        "betaC",
        dimensionSet(0, 0, 0, 0, 0),
        subDict(phase2Name_)
    ),
    TRef_
    (
        "TRef",
        dimensionSet(0, 0, 0, 1, 0),
        subDict(phase2Name_)
    ),
    CRef_
    (
        "CRef",
        dimensionSet(0, 0, 0, 0, 0),
        subDict(phase2Name_)
    ),

    U_(U)
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::solidificationSystem::read()
{
    const dictionary& phaseDict1_(subDict(phase1Name_));
    const dictionary& phaseDict2_(subDict(phase2Name_));

    rho1_ = phaseDict1_.lookupOrDefault<dimensionedScalar>("rho", dimensionedScalar("rho", dimDensity, 0.0));
    rho2_ = phaseDict2_.lookupOrDefault<dimensionedScalar>("rho", dimensionedScalar("rho", dimDensity, 0.0));

    Cp1_ = phaseDict1_.lookupOrDefault<dimensionedScalar>("Cp", dimensionedScalar("Cp", dimEnergy/dimMass/dimTemperature, 0.0));
    Cp2_ = phaseDict2_.lookupOrDefault<dimensionedScalar>("Cp", dimensionedScalar("Cp", dimEnergy/dimMass/dimTemperature, 0.0));

    kappa1_ = phaseDict1_.lookupOrDefault<dimensionedScalar>("kappa", dimensionedScalar("kappa", dimEnergy/dimLength/dimTemperature/dimTime, 0.0));
    kappa2_ = phaseDict2_.lookupOrDefault<dimensionedScalar>("kappa", dimensionedScalar("kappa", dimEnergy/dimLength/dimTemperature/dimTime, 0.0));

    mu1_ = phaseDict1_.lookupOrDefault<dimensionedScalar>("mu", dimensionedScalar("mu", dimMass/dimLength/dimTime, 0.0));
    mu2_ = phaseDict2_.lookupOrDefault<dimensionedScalar>("mu", dimensionedScalar("mu", dimMass/dimLength/dimTime, 0.0));

    D1_ = phaseDict1_.lookupOrDefault<dimensionedScalar>("D", dimensionedScalar("D", dimLength*dimLength/dimTime, 0.0));
    D2_ = phaseDict2_.lookupOrDefault<dimensionedScalar>("D", dimensionedScalar("D", dimLength*dimLength/dimTime, 0.0));

    DAS_ = phaseDict1_.lookupOrDefault<dimensionedScalar>("DAS", dimensionedScalar("DAS", dimLength, 0.0));

    betaT_ = phaseDict2_.lookupOrDefault<dimensionedScalar>("betaT", dimensionedScalar("betaT", dimless, 0.0));
    betaC_ = phaseDict2_.lookupOrDefault<dimensionedScalar>("betaC", dimensionedScalar("betaC", dimless, 0.0));

    TRef_ = phaseDict2_.lookupOrDefault<dimensionedScalar>("TRef", dimensionedScalar("TRef", dimTemperature, 0.0));
    CRef_ = phaseDict2_.lookupOrDefault<dimensionedScalar>("CRef", dimensionedScalar("CRef", dimless, 0.0));


    return true;
}


// ************************************************************************* //
