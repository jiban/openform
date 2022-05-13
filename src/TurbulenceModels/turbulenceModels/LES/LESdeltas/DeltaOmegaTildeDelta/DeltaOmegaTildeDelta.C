/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2022 Marian Fuchs
    Copyright (C) 2022 OpenCFD Ltd.
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

#include "DeltaOmegaTildeDelta.H"
#include "fvcCurl.H"
#include "maxDeltaxyz.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace LESModels
{
    defineTypeNameAndDebug(DeltaOmegaTildeDelta, 0);
    addToRunTimeSelectionTable(LESdelta, DeltaOmegaTildeDelta, dictionary);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::LESModels::DeltaOmegaTildeDelta::calcDelta()
{
    const fvMesh& mesh = turbulenceModel_.mesh();

    const volVectorField& U0 = turbulenceModel_.U();
    const volVectorField vorticity(fvc::curl(U0));
    const volVectorField nvecvort
    (
        vorticity
      / (
            max
            (
                mag(vorticity),
                dimensionedScalar("SMALL", dimless/dimTime, SMALL)
            )
        )
    );

    const cellList& cells = mesh.cells();
    const vectorField& cellC = mesh.cellCentres();
    scalarField hmax(cells.size());

    forAll(cells, celli)
    {
        scalar deltaMaxTmp = 0.0;

        const point& cc = cellC[celli];
        const labelList& cellPoints = mesh.cellPoints()[celli];
        const vector& nv = nvecvort[celli];

        // loop over all vertices
        for (const label pointi : cellPoints)
        {
            const point& pt = mesh.points()[pointi];
            deltaMaxTmp = max(mag(nv ^ (pt - cc)), deltaMaxTmp);
        }

        hmax[celli] = deltaCoeff_*Foam::sqrt(1.0/3.0)*2.0*deltaMaxTmp;
    }

    const label nD = mesh.nGeometricD();

    if (nD == 3)
    {
        delta_.primitiveFieldRef() = hmax;
    }
    else if (nD == 2)
    {
        WarningInFunction
            << "Case is 2D, LES is not strictly applicable\n"
            << endl;

        delta_.primitiveFieldRef() = hmax;
    }
    else
    {
        FatalErrorInFunction
            << "Case is not 3D or 2D, LES is not applicable"
            << exit(FatalError);
    }

    // Handle coupled boundaries
    delta_.correctBoundaryConditions();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::LESModels::DeltaOmegaTildeDelta::DeltaOmegaTildeDelta
(
    const word& name,
    const turbulenceModel& turbulence,
    const dictionary& dict
)
:
    LESdelta(name, turbulence),
    hmaxPtr_(nullptr),
    deltaCoeff_
    (
        dict.optionalSubDict(type() + "Coeffs").lookupOrDefault<scalar>
        (
            "deltaCoeff",
            1.025
        )
    ),
    requireUpdate_
    (
        dict.optionalSubDict(type() + "Coeffs").lookupOrDefault<bool>
        (
            "requireUpdate", true
        )
    )
{
    if (dict.optionalSubDict(type() + "Coeffs").found("hmax"))
    {
        // User-defined hmax
        hmaxPtr_ =
            LESdelta::New
            (
                IOobject::groupName("hmax", turbulence.U().group()),
                turbulence,
                dict.optionalSubDict(type() + "Coeffs"),
                "hmax"
            );
    }
    else
    {
        Info<< "Employing " << maxDeltaxyz::typeName << " for hmax" << endl;

        hmaxPtr_.reset
        (
            new maxDeltaxyz
            (
                IOobject::groupName("hmax", turbulence.U().group()),
                turbulence,
                dict.optionalSubDict(type() + "Coeffs")
            )
        );
    }

    calcDelta();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::LESModels::DeltaOmegaTildeDelta::read(const dictionary& dict)
{
    const dictionary& coeffsDict(dict.optionalSubDict(type() + "Coeffs"));

    coeffsDict.readIfPresent<scalar>("deltaCoeff", deltaCoeff_);
    coeffsDict.readIfPresent<bool>("requireUpdate", requireUpdate_);

    calcDelta();
}


void Foam::LESModels::DeltaOmegaTildeDelta::correct()
{
    if (turbulenceModel_.mesh().changing() && requireUpdate_)
    {
        hmaxPtr_->correct();
    }

    calcDelta();
}


// ************************************************************************* //
