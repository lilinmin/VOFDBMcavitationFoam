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
    Foam::bubbleCloud

Description
    A cloud of bubbles for the VoF-Lagrangian coupling.


SourceFiles
    bubbleCloudI.H
    bubbleCloud.C

\*---------------------------------------------------------------------------*/

#ifndef bubbleCloud_H
#define bubbleCloud_H

#include "Cloud.H"
#include "bubble.H"
#include "collisionModel.H"
#include "breakupModel.H"
#include "IOdictionary.H"
#include "CompactListList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class fvMesh;

/*---------------------------------------------------------------------------*\
                           Class bubbleCloud Declaration
\*---------------------------------------------------------------------------*/

class bubbleCloud
:
    public Cloud<bubble>
{
    // Private data

        const fvMesh& mesh_;

        const dimensionedVector& g_;
        
        IOdictionary dict_;
        IOdictionary particleProperties_;

        // bubble properties
        word phaseName_;
        
        scalar rhop_;
        scalar mup_;
        scalar sigma_;
        
        scalar p0_; //
        scalar pSat_; //added 7-6
       
        scalar gamma_;  //
        vector posP1_;
        scalar dP1_;
        vector UP1_;
        vector posP2_;
        scalar dP2_;
        vector UP2_;
        scalar tInjStart_;
        scalar tInjEnd_;
         
        //- Face zone IDs for postprocessing
        labelList faceZoneIDs_;

        //- Word list for face zone names
        wordList faceZoneNames_;

        //- Storage for output stream for postprocessing
        List<autoPtr<OFstream> > cloudDataFilesPtr_;
        
        //- Storage for bubble information for postprocessing
        List<List<DynamicList<scalar> > > dataDiameter_;
        List<List<DynamicList<scalar> > > dataNParticle_;
        List<List<DynamicList<vector> > > dataPosition_;
        
        // Auxiliary field for momentum source
        vectorField source_;
        vectorField Usource_;
        scalarField correctalpha1_;  // 7-99999
           
        // Momentum source for governing equation
        volVectorField momentumSource_;
        volVectorField USource_;
        volScalarField Addalpha1_;  // 7-99999

        //- Class for collision calculation
        collisionModel collision_;

        //- Class for breakup calculation
        breakupModel breakup_;

    // Private Member Functions

        //- Create file for bubble postprocessing
        void createbubbleDataFiles();

        //- Give a short summary of bubble cloud
        void info();

        //- Calculate mean diameter
        Foam::scalar Dij(const label i, const label j);

        //- Calculate maximum diameter
        Foam::scalar Dmax();


public:

    // Constructors

        //- Construct given mesh
        bubbleCloud
        (
            const fvMesh&,
            const dimensionedVector&, 
            const word& cloudName = "bubbleCloud",
            bool readFields = true
        );


    // Member Functions

        // Access

            inline const fvMesh& mesh() const;

            inline scalar rhop() const;
            inline scalar e() const;
            inline scalar mup() const;
            inline scalar sigma() const;
            
            inline scalar pSat() const;//added 7-6
             inline scalar RPEqn() const;// RP equation
            inline scalar p0() const;//added 7-6
            inline scalar gamma() const;//added 7-6
             
            // Field to store momentum source
            inline vectorField& source();
            inline vectorField& Usource();
            
            inline scalarField& correctalpha1();  //7-9
            // Return momentum source to carrier phase
            inline const volVectorField& momentumSource();
                 inline const volVectorField& USource(); 
            inline const volScalarField& Addalpha1();
             
            // Return faceZoneIDs for hitFace in bubble
            inline const labelList& faceZoneIDs();
            
            // Return processor data for bubble postprocessing
            inline List<List<DynamicList<scalar> > >& dataDiameter();
            inline List<List<DynamicList<vector> > >& dataPosition();
            inline List<List<DynamicList<scalar> > >& dataNParticle();
            
            

        // Edit

            //- Move the bubbles under the influence of the given
            //  gravitational acceleration
            void move();//bubbleCloud& cloud);
  
            //- Inject bubbles at a given position with a given
            //  velocity and diameter
            void inject
            (
                vector position,
                const scalar diameter,
                const vector velocity,
                bubbleCloud& cloud
            );
         
            void injector(bubbleCloud& cloud, bubble::trackingData &td); // For adding particles 
           
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "bubbleCloudI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
