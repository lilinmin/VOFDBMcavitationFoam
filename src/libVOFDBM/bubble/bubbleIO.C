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

#include "bubble.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const std::size_t Foam::bubble::sizeofFields
(
    sizeof(bubble) - sizeof(particle)
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::bubble::bubble
(
    const polyMesh& mesh,
    Istream& is,
    bool readFields,
    bool newFormat
)
:
    particle(mesh, is, readFields, newFormat)
{
    if (readFields)
    {
        if (is.format() == IOstream::ASCII)
        {
            d_ = readScalar(is);
            is >> U_>> Rdt_ >> R_ >> R0_ >> pB_ >> pG_ >> pG0_ >> F0Old_ >> F1Old_;
            nParticle_ = readScalar(is);
            y_ = readScalar(is);
            yDot_ = readScalar(is);
        }
        else
        {
            is.read(reinterpret_cast<char*>(&d_), 
            sizeof(d_)
+ sizeof(U_)
+ sizeof(Rdt_)
+ sizeof(R_)
+ sizeof(R0_)
+ sizeof(pB_)
+ sizeof(pG_)
+ sizeof(pG0_)
+ sizeof(F0Old_)
+ sizeof(F1Old_)
+ sizeof(nParticle_)
+ sizeof(y_)
+ sizeof(yDot_));
              
            //sizeofFields); //modified try
        }
    }

    is.check(FUNCTION_NAME);
}


void Foam::bubble::readFields(Cloud<bubble>& c)
{
    bool valid = c.size();

    particle::readFields(c);

    IOField<scalar> d(c.fieldIOobject("d", IOobject::MUST_READ), valid);
    c.checkFieldIOobject(c, d);

    IOField<vector> U(c.fieldIOobject("U", IOobject::MUST_READ), valid);
    c.checkFieldIOobject(c, U);

IOField<scalar> Rdt(c.fieldIOobject("Rdt", IOobject::MUST_READ), valid);
c.checkFieldIOobject(c, Rdt);
IOField<scalar> R(c.fieldIOobject("R", IOobject::MUST_READ), valid);
c.checkFieldIOobject(c, R);
IOField<scalar> R0(c.fieldIOobject("R0", IOobject::MUST_READ), valid);
c.checkFieldIOobject(c, R0);
IOField<scalar> pB(c.fieldIOobject("pB", IOobject::MUST_READ), valid);
c.checkFieldIOobject(c, pB);
IOField<scalar> pG(c.fieldIOobject("pG", IOobject::MUST_READ), valid);
c.checkFieldIOobject(c, pG);
IOField<scalar> pG0(c.fieldIOobject("pG0", IOobject::MUST_READ), valid);
c.checkFieldIOobject(c, pG0);
IOField<scalar> F0Old(c.fieldIOobject("F0Old", IOobject::MUST_READ), valid);
c.checkFieldIOobject(c, F0Old);
IOField<scalar> F1Old(c.fieldIOobject("F1Old", IOobject::MUST_READ), valid);
c.checkFieldIOobject(c, F1Old);

    IOField<scalar> nParticle(c.fieldIOobject("nParticle", IOobject::MUST_READ), valid);
    c.checkFieldIOobject(c, nParticle);

    IOField<scalar> y(c.fieldIOobject("y", IOobject::MUST_READ), valid);
    c.checkFieldIOobject(c, y);

    IOField<scalar> yDot(c.fieldIOobject("yDot", IOobject::MUST_READ), valid);
    c.checkFieldIOobject(c, yDot);

    label i = 0;
    forAllIter(Cloud<bubble>, c, iter)
    {
        bubble& p = iter();

        p.d_ = d[i];
        p.U_ = U[i];
        
p.Rdt_ = Rdt[i];
p.R_ = R[i];
p.R0_ = R0[i];
p.pB_ = pB[i];
p.pG_ = pG[i];
p.pG0_ = pG0[i];
p.F0Old_ = F0Old[i];
p.F1Old_ = F1Old[i];  

        p.nParticle_ = nParticle[i];
        p.y_ = y[i];
        p.yDot_ = yDot[i];
        i++;
    }
}


void Foam::bubble::writeFields(const Cloud<bubble>& c)
{
    particle::writeFields(c);

    label np = c.size();

    IOField<scalar> d(c.fieldIOobject("d", IOobject::NO_READ), np);
    IOField<vector> U(c.fieldIOobject("U", IOobject::NO_READ), np);
    
IOField<scalar> Rdt(c.fieldIOobject("Rdt", IOobject::NO_READ), np);
IOField<scalar> R(c.fieldIOobject("R", IOobject::NO_READ), np);
IOField<scalar> R0(c.fieldIOobject("R0", IOobject::NO_READ), np);
IOField<scalar> pB(c.fieldIOobject("pB", IOobject::NO_READ), np);
IOField<scalar> pG(c.fieldIOobject("pG", IOobject::NO_READ), np);
IOField<scalar> pG0(c.fieldIOobject("pG0", IOobject::NO_READ), np);
IOField<scalar> F0Old(c.fieldIOobject("F0Old", IOobject::NO_READ), np);
IOField<scalar> F1Old(c.fieldIOobject("F1Old", IOobject::NO_READ), np);
     
    IOField<scalar> nParticle(c.fieldIOobject("nParticle", IOobject::NO_READ), np);
    IOField<scalar> y(c.fieldIOobject("y", IOobject::NO_READ), np);
    IOField<scalar> yDot(c.fieldIOobject("yDot", IOobject::NO_READ), np);

    label i = 0;
    forAllConstIter(Cloud<bubble>, c, iter)
    {
        const bubble& p = iter();

        d[i] = p.d_;
        U[i] = p.U_;
        
Rdt[i] = p.Rdt_;
R[i] = p.R_;
R0[i] = p.R0_;
pB[i] = p.pB_;
pG[i] = p.pG_;
pG0[i] = p.pG0_;
F0Old[i] = p.F0Old_;
F1Old[i] =p.F1Old_;
        
        nParticle[i] = p.nParticle_;
        y[i] = p.y_;
        yDot[i] = p.yDot_;
        i++;
    }

    d.write(np > 0);
    U.write(np > 0);
    
Rdt.write(np > 0);
R.write(np > 0);
R0.write(np > 0);
pB.write(np > 0);
pG.write(np > 0);
pG0.write(np > 0);
F0Old.write(np > 0);
F1Old.write(np > 0);
    
    nParticle.write(np > 0);
    y.write(np > 0);
    yDot.write(np > 0);
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const bubble& p)
{
    if (os.format() == IOstream::ASCII)
    {
        os  << static_cast<const particle&>(p)
            << token::SPACE << p.d_
            << token::SPACE << p.U_
            
            << token::SPACE << p.Rdt_
            << token::SPACE << p.R_
            << token::SPACE << p.R0_
            << token::SPACE << p.pB_
            << token::SPACE << p.pG_
            << token::SPACE << p.pG0_
            << token::SPACE << p.F0Old_
            << token::SPACE << p.F1Old_
            
            << token::SPACE << p.nParticle_;
    }
    else
    {
        os  << static_cast<const particle&>(p);
        os.write
        (
            reinterpret_cast<const char*>(&p.d_),
sizeof(p.d_)
+ sizeof(p.U_)
+ sizeof(p.Rdt_)
+ sizeof(p.R_)
+ sizeof(p.R0_)
+ sizeof(p.pB_)
+ sizeof(p.pG_)
+ sizeof(p.pG0_)
+ sizeof(p.F0Old_)
+ sizeof(p.F1Old_)
+ sizeof(p.nParticle_)
+ sizeof(p.y_)
+ sizeof(p.yDot_)
 
    //bubble::sizeofFields
        );
    }

    os.check(FUNCTION_NAME);
    return os;
}


// ************************************************************************* //
