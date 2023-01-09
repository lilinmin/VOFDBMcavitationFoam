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

Class
    Foam::bubble

Description
    For solving the RP equation.
  
SourceFiles
    bubbleI.H
    bubble.C
    bubbleIO.C


\*---------------------------------------------------------------------------*/

//Solve R-P equation for bubbles
scalar hrk = 1.0E-10;
scalar h = 0.5*hrk;
scalar dD = 1.0/(2.0+sqrt(2.0));
scalar time_rk = 0;
scalar sigma = cloud.sigma();
//scalar pV = cloud.pSat();
scalar substeps = 0;
 
while (time_rk < trackTime)
{ 
 if (time_rk+h > trackTime)
 {
    h = trackTime-time_rk;
 }
// f
scalar F0 = Rdt_;
scalar F1 = -1.5*(Rdt_*Rdt_/R_)+((pB_-pc)/(rhoc*R_))-4.0*(nuc*rhoc)*(Rdt_/(R_*R_))-((2.0*sigma)/(rhoc*R_*R_));
// f0
scalar F00 = F0;
scalar F01 = F1;
// the jacobian dF/dX
scalar J00 = 0.0;
scalar J01 = 1.0;
scalar J10 = 1.5*((Rdt_*Rdt_)/(R_*R_))-((pB_-pc)/(rhoc*R_*R_))+
((8.0*(nuc*rhoc)*Rdt_)/(R_*R_*R_))+((4.0*sigma)/(rhoc*R_*R_*R_));
scalar J11 = -3.0*Rdt_/R_-4.0*((nuc*rhoc)/(R_*R_));
// W = I - h dD J
scalar W00 = 1.0 - h*dD*J00;
scalar W01 = 0.0 - h*dD*J01;
scalar W10 = 0.0 - h*dD*J10;
scalar W11 = 1.0 - h*dD*J11;
// inv W
scalar ff = 1.0/(W00*W11-W01*W10);
scalar invW00 = ff*W11;
scalar invW01 = -ff*W01;
scalar invW10 = -ff*W10;
scalar invW11 = ff*W00;
// the time derivative of function, dF/dt
// estimated numerically
scalar T0 = (F0-F0Old_)/h;
scalar T1 = (F1-F1Old_)/h;
// k1
scalar tmp0 = F00+h*dD*T0;
scalar tmp1 = F01+h*dD*T1;
scalar k10 = invW00*tmp0+invW01*tmp1;
scalar k11 = invW10*tmp0+invW11*tmp1;
// y_n + 0.5 h k1
scalar R1 = R_+0.5*h*k10;
scalar Rdt1 = Rdt_+0.5*h*k11;
// f1
scalar F10 = Rdt1;
scalar F11 = -1.5*(Rdt1*Rdt1/R1)+((pB_ -
pc)/(rhoc*R1))-4.0*(nuc*rhoc)*(Rdt1/(R1*R1))-((2.0*sigma)/(rhoc*R1*R1));
// k2
tmp0 = F10-k10;
tmp1 = F11-k11;
scalar k20 = invW00*tmp0+invW01*tmp1+k10;
scalar k21 = invW10*tmp0+invW11*tmp1+k11;
// y_n+1
scalar R2 = R_+h*k20;
scalar Rdt2 = Rdt_+h*k21;
// f2
scalar F20 = Rdt2;
scalar F21 = -1.5*(Rdt2*Rdt2/R2)+((pB_ -
pc)/(rhoc*R2))-4.0*(nuc*rhoc)*(Rdt2/(R2*R2))-((2.0*sigma)/(rhoc*R2*R2));

// k3
scalar b = 6.0 + sqrt(2.0);
tmp0 = F20-b*(k20-F10)-2.0*(k10-F00)+h*dD*T0;
tmp1 = F21-b*(k21-F11)-2.0*(k11-F01)+h*dD*T1;
scalar k30 = invW00*tmp0+invW01*tmp1;
scalar k31 = invW10*tmp0+invW11*tmp1;
// error
scalar err0 = fabs((h/6.0)*(k10-2.0*k20+k30));
scalar err1 = fabs((h/6.0)*(k11-2.0*k21+k31));
scalar err = max(err0,err1);

time_rk += h;
substeps += 1;
//scalar absTol = 1.0E-5;

if (err > 1.0E-08)
{
h = min(0.5*hrk, 0.5*h);
}
else
{
h = max(1.5*hrk, 1.5*h);
}
// update bubble state
F0Old_ = F0;
F1Old_ = F1;
R_ = R2;
d_ = 2.0*R_;
Rdt_ = Rdt2;

if (R_ < R0_)
{
  pG_ = pG0_*pow(R0_/R_,3.0);
}
else
{
  pG_ = pG0_*pow(R0_/R_,3.0*1.4);
}
  pB_ = pG_ + pV;

}
// End of R-P equation 

