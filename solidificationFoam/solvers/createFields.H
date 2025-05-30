Info<< "Reading field p_rgh\n" << endl;
volScalarField p_rgh
(
    IOobject
    (
        "p_rgh",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

solidificationSystem fluid(U);

Info<< "Reading field T\n" << endl;
volScalarField T
(
    IOobject
    (
        "T",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field C\n" << endl;
volScalarField C
(
    IOobject
    (
        "C",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

#include "createAuxilaryFields.H"

// Mixure density for shrinkage force
volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh
    ),
    fluid.rho()
);

#include "compressibleCreatePhi.H"

// Kinematic density for buoyancy force
volScalarField rhok
(
    IOobject
    (
        "rhok",
        runTime.timeName(),
        mesh
    ),
    fluid.rho2()
   *(1.0 - fluid.betaT()*(T - fluid.TRef()) - fluid.betaC()*(C2 - fluid.CRef()))
);


#include "readGravitationalAcceleration.H"
#include "readhRef.H"
#include "gh.H"

volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    p_rgh + rhok*gh
);

label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell
(
    p,
    p_rgh,
    pimple.dict(),
    pRefCell,
    pRefValue
);

if (p_rgh.needReference())
{
    p += dimensionedScalar
    (
        "p",
        p.dimensions(),
        pRefValue - getRefCellValue(p, pRefCell)
    );
}

mesh.setFluxRequired(p_rgh.name());

#include "createMRF.H"
#include "createFvOptions.H"

{
    Info<< "Detected cellZones:" << endl;
    const cellZoneMesh& cz = mesh.cellZones();
    const List<word>& zoneNames = cz.cellZoneNames();
    forAll(zoneNames, zi)
    {
        const word& name = zoneNames[zi];
        label nCells = cz[name].size();
        Info<< "  - " << name << " (" << nCells << " cells)" << endl;
    }
}
{
    Info<< "Detected boundary patches:" << endl;
    const fvPatchList& pList = mesh.boundaryMesh();
    const List<word>& patchNames = pList.patchNames();
    forAll(patchNames, pi)
    {
        Info<< "  - " << patchNames[pi] << endl;
    }
}
if (mesh.cellZones().found("stator"))
{
    labelList const& cells = mesh.cellZones()["stator"];
    Info<< "Creating mask for zone stator with " << cells.size() << " cells" << endl;
}
Info<< "Finished dynamic zone/patch detection\n" << endl;
