/*

 Copyright 2010 Brain Research Institute/National ICT Australia (NICTA), Melbourne, Australia

 Modified by Thomas G Close, 10/4/11 from KMlocal.cpp
 written by David Mount from the University of Maryland, 05/14/04
 as part of the Kmeans package available at
 'http://www.cs.umd.edu/~mount/Projects/KMeans/'
 under the GNU General public licence.

 This file is part of Bayesian Tractlet Sampling (BTS).

 BTS is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 BTS is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with BTS.  If not, see <http://www.gnu.org/licenses/>.

 */

#include "bts/math/kmeans.h"				// KMlocal includes
namespace BTS {
    
    namespace Math {
        
        //----------------------------------------------------------------------
        //  execute - execute the clustering algorithm
        //  	This function executes the clustering algorithm.  See the file
        //  	KMlocal.h for a description of this algorithm.
        //----------------------------------------------------------------------
        
        ::KMfilterCenters KMeansBase::execute()		// execute the algorithm
        {
            reset();					// resets everythinnessg
            while (!isDone()) {				// while not done
                beginRun();				// start a new run
                do {					// do while run is not done
                    beginStage();			// start of stage processing
                    KMalg method = selectMethod();    // select a method
                    switch (method) {			// apply one stage
                        case LLOYD:				// Lloyd's algorithm
                            curr.lloyd1Stage();
                            break;
                        case SWAP:				// swap heuristic
                            curr.swap1Stage();
                            break;
                        case RANDOM:			// get random centers
                            curr.genRandom();
                            break;
                        default:				// shouldn't come here
                            assert(false);
                            break;
                    }
                    endStage();				// end of stage processing
                } while (!isRunDone());			// while run is not done
                endRun();				// end of run processing
                tryAcceptance();			// accept if appropriate
            }
            return best;				// return best solution
        }
    
    }
}
