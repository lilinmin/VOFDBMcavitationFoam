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

#include "bubbleCloud.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTemplateTypeNameAndDebug(Cloud<bubble>, 0);
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::bubble::move
(
    bubbleCloud& cloud,
    trackingData& td, 
    const scalar trackTime
)
{
    td.switchProcessor = false;
    td.keepParticle = true;
     

    while (td.keepParticle && !td.switchProcessor && stepFraction() < 1)
    {
        if (debug)
        {
            Info<< "Time = " << mesh().time().timeName()
                << " trackTime = " << trackTime
                << " steptFraction() = " << stepFraction() << endl;
        }

        const scalar sfrac = stepFraction();

        const scalar f = 1 - stepFraction();
        trackToAndHitFace(f*trackTime*U_, f, cloud, td);

        const scalar dt = (stepFraction() - sfrac)*trackTime;

        const tetIndices tetIs = this->currentTetIndices();
        scalar rhoc = td.rhoInterp().interpolate(this->coordinates(), tetIs);
        vector Uc = td.UInterp().interpolate(this->coordinates(), tetIs);
        scalar nuc = td.nuInterp().interpolate(this->coordinates(), tetIs);
        //scalar alpha1P = td.alpha1Interp().interpolate(this->coordinates(), tetIs);//
        scalar pc = td.pInterp().interpolate(this->coordinates(), tetIs);//

        scalar rhop = cloud.rhop();
        scalar magUr = mag(Uc - U_);

        //    label cellI = cell();

        scalar ReFunc = 1.0;
        scalar Re = magUr*d_/nuc;
        
        if (Re >= 0.1 && Re <= 1000.0)
        {
            ReFunc += 0.15*pow(Re, 0.687);
        }
        if (Re > 1000.0)
        { 
            ReFunc = Re/24.0*0.44;
        }

        scalar Dc = (24.0*nuc/d_)*ReFunc*(3.0/4.0)*(rhoc/(d_*rhop));
        
        // calculate momentum
        scalar m = rhop*(4.0/3.0)*constant::mathematical::pi*pow(d_/2.0, 3.0);
        vector i1 = U_*m;
        
        U_ = (U_ + dt*(Dc*Uc + (1.0 - rhoc/rhop)*td.g()))/(1.0 + dt*Dc);

//For bubble size variation
    scalar pV = cloud.pSat();
    
    if (cloud.RPEqn() == 1)
    {
       #include "RPEqn.H"
    }
    else 
    {
       if( pc > 1.25*pV )
       {
          d_ = max(0.000001, d_- 0.01*sqrt(0.66*(pc - 1.25*pV)/rhoc)*dt);
       }
       else
       {
          d_ = min(0.02, d_+ 0.5*sqrt(0.66*(1.25*pV - pc)/rhoc)*dt); 
       } 
    }
    
//LPT to VOF by Linmin Li, it is not suggested for phasechange cases
bool nearInterface = false;

vector posP = position(); 
label cellP = mesh().findCell(posP);
        if (mag(i1) > 0.0)
        { 
              cloud.source()[cellP] += - (U_*m - i1 - m * (1.0 - rhoc/rhop) * td.g() * dt); 
        }
scalar alpha1P = td.alpha1Interp().interpolate(posP, cellP);
scalar distanceP = d_/2.0;
for (int j=0; j<6; j++)
{
  int sign=pow(-1,j);
  int coord=j/2;
  vector possurfP=posP;
  possurfP.component(coord) += sign * distanceP ;
  label cellsurfP = mesh().findCell(possurfP);
  if (cellsurfP != (-1))
   {
    scalar alpha1surfP = td.alpha1Interp().interpolate(possurfP, cellsurfP);
    if (alpha1surfP >= 0.5)
     {
       nearInterface = true;
     }
   }
}

scalar VPleft;
scalar VdispoCellP;
 
if (nearInterface)
{
   VPleft=pow(d_,3)*Foam::constant::mathematical::pi/6.0;
   VdispoCellP = (1.0-alpha1P)*mesh().cellVolumes()[cellP];

   if (VPleft <= VdispoCellP)
   //fill only the eulerian cell at the particle position
   {
       cloud.correctalpha1()[cellP] += VPleft/mesh().cellVolumes()[cellP]; 
       cloud.Usource()[cellP] = (U_-Uc)/dt; 
   }
   else
  // start to fill the neighbour cells
   {
      cloud.correctalpha1()[cellP] = 1.0; 
      
      cloud.Usource()[cellP] = (U_-Uc)/dt; 
      
      // Create dynamic list to store all neighbours
      DynamicList<label> cells(0);
      // Get neighbours of particle Cell
      cells.append(mesh().cellCells()[cellP]);
 
        // Loop over all current neighbours
        forAll(cells, i)
        {
              label cellIn = cells[i]; 
            
         }
   }
   td.keepParticle = false;
}  
//End of LPT to VOF transition 


if (d_<=0.000001 || d_>=0.001 || alpha1P>0.8)
{
     td.keepParticle = false;
}
        // Check if particle passed given faceZones
        if (face() > -1)
        {
            const faceZoneMesh& fzm = mesh().faceZones();
            
            const labelList fzIDs = cloud.faceZoneIDs();
            
            forAll(fzIDs, i)
            {
                const faceZone& fz = fzm[fzIDs[i]];

                forAll(fz, j)
                {                
                    if (fz[j] == face())
                    {
                        cloud.dataDiameter()[Pstream::myProcNo()][i].append(d_);
                        cloud.dataPosition()[Pstream::myProcNo()][i].append(position());
                        cloud.dataNParticle()[Pstream::myProcNo()][i].append(nParticle_);
                    }
                }
            }
        }
         
    }

    return td.keepParticle;
}

/* Rayleigh-Plesset: initial equilibrium radius */
//
bool Foam::bubble::initialEquilibriumRadius
(
   bubbleCloud& cloud,
   trackingData& td
)
{
  td.switchProcessor = false;
  td.keepParticle = true;

// remember which cell the parcel is in
// since this will change if a face is hit
    const tetIndices tetIs = this->currentTetIndices();
    scalar pc = td.pInterp().interpolate(this->coordinates(), tetIs);//

/// initial equlibirium radius
scalar p0 = cloud.p0();
scalar pV = cloud.pSat();
scalar sigma = cloud.sigma();
scalar gamma = cloud.gamma();
scalar R = 0.5*d_;
R_ = 0.5*d_;
Rdt_ = 0.0;
F0Old_ = 0.0;
F1Old_ = 0.0;
// initial guess R0 equal to R
scalar R0 = R;
scalar R0new = R;
scalar R0bar = R;
scalar error = 1;

scalar a = (2.0*sigma)/(p0-pV);
int maxIter = 2000;
int iter = 0;
while ((error > 1.0E-06) && (iter < maxIter))
{
if (R < R0)
{
  gamma = 1.4;
}
else
{
  gamma = 1.0;
}
  scalar b = (pow(R,3.0*gamma)/(p0-pV))*(pV-((2.0*sigma)/R)-pc);
  scalar f = pow(R0,3.0*gamma) + a*pow(R0,3.0*gamma-1.0) + b;
  scalar df = 3.0*gamma*pow(R0,3*gamma-1.0) +
  a*(3.0*gamma-1.0)*pow(R0,3.0*gamma-2.0);
  R0bar = R0 - (f/df);
if (R < R0bar)
{
  gamma = 1.4;
}
else
{
  gamma = 1.0;
}
  b = (pow(R,3.0*gamma)/(p0-pV))*(pV-((2.0*sigma)/R)-pc);
  scalar fbar = pow(R0bar,3.0*gamma) + a*pow(R0bar,3.0*gamma-1.0) + b;
  R0new = R0 - (f*f)/(df*(f-fbar));
  R0new = 0.9*R0 + 0.1*R0new;
  error = fabs((R0new-R0)/R0);
  R0 = R0new;
  iter++;
}
R0_ = R0;
pG0_ = p0 + 2.0*sigma/R0_ - pV;
if (R < R0)
{
  gamma = 1.4;
}
else
{
  gamma = 1.0;
}
  pG_ = pG0_*pow(R0_/R_,3.0);
  pB_ = pG_ + pV;

  return true;
}
//added7-6

bool Foam::bubble::hitPatch(bubbleCloud&, trackingData&)
{
    return false;
}


void Foam::bubble::hitProcessorPatch
(
    bubbleCloud&,
    trackingData& td
)
{
    td.switchProcessor = true;
}


void Foam::bubble::hitWallPatch(bubbleCloud& cloud, trackingData&)
{
    const vector nw = normal();

    scalar Un = U_ & nw;
    vector Ut = U_ - Un*nw;

    if (Un > 0)
    {
        U_ -= (1.0 + 0.8)*Un*nw;
    }

    U_ -= 0.09*Ut;
}


void Foam::bubble::transformProperties(const tensor& T)
{
    particle::transformProperties(T);
    U_ = transform(T, U_);
}


void Foam::bubble::transformProperties(const vector& separation)
{
    particle::transformProperties(separation);
}


// ************************************************************************* //
