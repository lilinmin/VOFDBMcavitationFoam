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
    Spherical bubble class with optional two-way coupling with the
    continuous phase, bubble-bubble collision, and breakup..

SourceFiles
    bubbleI.H
    bubble.C
    bubbleIO.C

\*---------------------------------------------------------------------------*/

#ifndef bubble_H
#define bubble_H

#include "particle.H"
#include "IOstream.H"
#include "autoPtr.H"
#include "interpolationCellPoint.H"
#include "interpolationCell.H" //
#include "contiguous.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class bubbleCloud;


// Forward declaration of friend functions and operators

class bubble;

Ostream& operator<<(Ostream&, const bubble&);


/*---------------------------------------------------------------------------*\
                           Class bubble Declaration
\*---------------------------------------------------------------------------*/

class bubble
:
    public particle
{
    // Private data

        //- Diameter
        scalar d_;

        //- Velocity of particle
        vector U_;

//- temporal derivative of the Radius
scalar Rdt_;
//- Radius
scalar R_;
//- Initial equilibirium radius
scalar R0_;
//- bubble inside pressure
scalar pB_;
//- dissolved gas pressure
scalar pG_;
scalar pG0_;
scalar F0Old_;
scalar F1Old_;
//added
        //- Number of particles in parcel
        scalar nParticle_;

        //- Spherical deviation
	scalar y_;

	//- Rate of change of spherical deviation
	scalar yDot_;

public:

    friend class Cloud<bubble>;

    //- Class used to pass tracking data to the trackToFace function
    class trackingData
    :
        public particle::trackingData
    {
        // Interpolators for continuous phase fields

            const interpolationCellPoint<scalar>& rhoInterp_;
            const interpolationCellPoint<vector>& UInterp_;
            const interpolationCellPoint<scalar>& nuInterp_;
            const interpolationCell<scalar>& alpha1Interp_; //
            const interpolationCellPoint<scalar>& pInterp_; // 7-6

        //- Local gravitational or other body-force acceleration
        const vector& g_;


    public:

        // Constructors

            inline trackingData
            (
                const bubbleCloud& spc,
                const interpolationCellPoint<scalar>& rhoInterp,
                const interpolationCellPoint<vector>& UInterp,
                const interpolationCellPoint<scalar>& nuInterp,
                const interpolationCell<scalar>& alpha1Interp, //     
                const interpolationCellPoint<scalar>& pInterp_, // 7-6
                const vector& g
            );


        // Member functions

            inline const interpolationCellPoint<scalar>& rhoInterp() const;

            inline const interpolationCellPoint<vector>& UInterp() const;

            inline const interpolationCellPoint<scalar>& nuInterp() const;
            
            inline const interpolationCell<scalar>& alpha1Interp() const; //    
            inline const interpolationCellPoint<scalar>& pInterp() const; // 7-6

            inline const vector& g() const;
    };


    // Static data members

        //- Size in bytes of the fields
        static const std::size_t sizeofFields;


    // Constructors

        //- Construct from components
        inline bubble
        (
            const polyMesh& mesh,
            const vector& position,
            const label celli,
            const scalar d,
            const vector& U
        );

        //- Construct from Istream
        bubble
        (
            const polyMesh& mesh,
            Istream& is,
            bool readFields = true,
            bool newFormat = true
        );

        //- Construct and return a clone
        virtual autoPtr<particle> clone() const
        {
            return autoPtr<particle>(new bubble(*this));
        }

        //- Factory class to read-construct particles used for
        //  parallel transfer
        class iNew
        {
            const polyMesh& mesh_;

        public:

            iNew(const polyMesh& mesh)
            :
                mesh_(mesh)
            {}

            autoPtr<bubble> operator()(Istream& is) const
            {
                return autoPtr<bubble>
                (
                    new bubble(mesh_, is, true)
                );
            }
        };


    // Member Functions

        // Access

            //- Return diameter
            inline scalar& d();

            //- Return velocity
            inline vector& U();

            //- Return number of particles in parcel
            inline scalar& nParticle();

            //- Return spherical deviation
            inline scalar& y();

            //- Return change of rate of spherical deviation
            inline scalar& yDot();


        // Tracking
bool initialEquilibriumRadius(bubbleCloud& cloud, trackingData& td); //add 7-6

            //- Move
            bool move(bubbleCloud&, trackingData&, const scalar);


        // Patch interactions

            //- Overridable function to handle the particle hitting a patch
            //  Executed before other patch-hitting functions
            bool hitPatch(bubbleCloud& cloud, trackingData& td);

            //- Overridable function to handle the particle hitting a
            //  processorPatch
            void hitProcessorPatch(bubbleCloud& cloud, trackingData& td);

            //- Overridable function to handle the particle hitting a wallPatch
            void hitWallPatch(bubbleCloud& cloud, trackingData& td);

            //- Transform the physical properties of the particle
            //  according to the given transformation tensor
            virtual void transformProperties(const tensor& T);

            //- Transform the physical properties of the particle
            //  according to the given separation vector
            virtual void transformProperties(const vector& separation);


    // I-O

        static void readFields(Cloud<bubble>& c);

        static void writeFields(const Cloud<bubble>& c);


    // Ostream Operator

        friend Ostream& operator<<(Ostream&, const bubble&);
};

//- Contiguous data for solidParticle
template<> struct is_contiguous<bubble> : std::true_type {};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "bubbleI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
