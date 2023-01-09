/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::bubble::trackingData::trackingData
(
    const bubbleCloud& spc,
    const interpolationCellPoint<scalar>& rhoInterp,
    const interpolationCellPoint<vector>& UInterp,
    const interpolationCellPoint<scalar>& nuInterp,
    const interpolationCell<scalar>& alpha1Interp, //
    const interpolationCellPoint<scalar>& pInterp, //
    const vector& g
)
:
    particle::trackingData(spc),
    rhoInterp_(rhoInterp),
    UInterp_(UInterp),
    nuInterp_(nuInterp),
    alpha1Interp_(alpha1Interp), //
    pInterp_(pInterp),//
    g_(g)
{}


inline Foam::bubble::bubble
(
    const polyMesh& mesh,
    const vector& position,
    const label celli,
    const scalar d,
    const vector& U
)
:
    particle(mesh, position, celli),
    d_(d),
    U_(U),
    nParticle_(1),
    y_(0.0),
    yDot_(0.0)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline const Foam::interpolationCellPoint<Foam::scalar>&
Foam::bubble::trackingData::rhoInterp() const
{
    return rhoInterp_;
}


inline const Foam::interpolationCellPoint<Foam::vector>&
Foam::bubble::trackingData::UInterp() const
{
    return UInterp_;
}


inline const Foam::interpolationCellPoint<Foam::scalar>&
Foam::bubble::trackingData::nuInterp() const
{
    return nuInterp_;
}

//added
inline const Foam::interpolationCell<Foam::scalar>&
Foam::bubble::trackingData::alpha1Interp() const
{
    return alpha1Interp_;
}
//end
//added
inline const Foam::interpolationCellPoint<Foam::scalar>&
Foam::bubble::trackingData::pInterp() const
{
    return pInterp_;
}
//end
inline const Foam::vector& Foam::bubble::trackingData::g() const
{
    return g_;
}


inline Foam::scalar& Foam::bubble::d()
{
    return d_;
}


inline Foam::vector& Foam::bubble::U()
{
    return U_;
}


inline Foam::scalar& Foam::bubble::nParticle()
{
    return nParticle_;
}


inline Foam::scalar& Foam::bubble::y()
{
    return y_;
}


inline Foam::scalar& Foam::bubble::yDot()
{
    return yDot_;
}


// ************************************************************************* //
