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

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline const Foam::fvMesh& Foam::bubbleCloud::mesh() const
{
    return mesh_;
}

inline Foam::scalar Foam::bubbleCloud::rhop() const
{
    return rhop_;
}

inline Foam::scalar Foam::bubbleCloud::mup() const
{
    return mup_;
}

inline Foam::scalar Foam::bubbleCloud::sigma() const
{
    return sigma_;
}

//added
inline Foam::scalar Foam::bubbleCloud::p0() const
{
  return p0_;
}

inline Foam::scalar Foam::bubbleCloud::pSat() const
{
    return pSat_;
}
inline Foam::scalar Foam::bubbleCloud::RPEqn() const
{
    return RPEqn_;
}

inline Foam::scalar Foam::bubbleCloud::gamma() const
{
    return gamma_;
}

//7-6
inline Foam::vectorField& Foam::bubbleCloud::source()
{
    return source_;
}
inline Foam::vectorField& Foam::bubbleCloud::Usource()
{
    return Usource_;
}
inline Foam::scalarField& Foam::bubbleCloud::correctalpha1()
{
     return correctalpha1_;
}
inline const Foam::volVectorField& Foam::bubbleCloud::momentumSource()
{
    return momentumSource_;
}
inline const Foam::volVectorField& Foam::bubbleCloud::USource()
{
    return USource_;
}
inline const Foam::volScalarField& Foam::bubbleCloud::Addalpha1()
{ 
    return Addalpha1_;
}

inline const Foam::labelList& Foam::bubbleCloud::faceZoneIDs()
{
    return faceZoneIDs_;
}

inline Foam::List<Foam::List<Foam::DynamicList<Foam::scalar> > >& Foam::bubbleCloud::dataDiameter()
{
    return dataDiameter_;
}

inline Foam::List<Foam::List<Foam::DynamicList<Foam::scalar> > >& Foam::bubbleCloud::dataNParticle()
{
    return dataNParticle_;
}

inline Foam::List<Foam::List<Foam::DynamicList<Foam::vector> > >& Foam::bubbleCloud::dataPosition()
{
    return dataPosition_;
}
//added

//end
// ************************************************************************* //
