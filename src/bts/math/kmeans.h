/*

    Copyright 2010 Brain Research Institute/National ICT Australia (NICTA), Melbourne, Australia

    Modified by Thomas G Close, 10/4/11 from KMlocal.h
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




#ifndef __bts_math_kmeans_h__
#define __bts_math_kmeans_h__

//----------------------------------------------------------------------
//  basic includes
//----------------------------------------------------------------------
#include <cmath>				// math includes (exp, log)
#include "k_means/KMeans.h"				// kmeans includes
#include "k_means/KMdata.h"				// data points
#include "k_means/KMfilterCenters.h"			// centers
#include "k_means/KMterm.h"				// termination conditions
#include "k_means/KMrand.h"				// random number generation
#include "k_means/KMlocal.h"


namespace BTS {

  namespace Math {

    //------------------------------------------------------------------------
    //  Generic Local Search Overview
    //  KMlocal provides a generic algorithm for k-means clustering by
    //  local search.  This generic algorithm is controlled by a class,
    //  KMlocal, whose methods provide the particular elements of each
    //  algorithm.
    //
    //  Overview:
    //  ---------
    //  The generic algorithm begins by generating an initial solution
    //  "curr" and saving it in "best".  These objects are local to the
    //  KMlocal structure.   The value of "curr" reflects the current
    //  solution and "best" reflects the best solution seen so far.
    //
    //  The algorithm consists of some number of basic iterations,
    //  called "stages".  Each stage involves the execution of one step
    //  of either the swap heuristic or Lloyd's algorithm.
    //
    //  Stages are grouped into "runs".  Intuitively, a run involves a
    //  (small) number of stages in search of a better solution.  A run
    //  might end, say, because a better solution was found or a fixed
    //  number of stages have been performed without any improvement.
    //
    //  After a run is finished, we check to see whether we want to
    //  "accept" the solution.  Presumably this happens if the cost is
    //  lower, but it may happen even if the cost is inferior in other
    //  circumstances (e.g., as part of a simulated annealing approach).
    //  Accepting a solution means copying the current solution to the
    //  saved solution.  In some cases, the acceptance may involve
    //  reseting the current solution to a random solution.
    //
    //  Generic Algorithm:
    //  ------------------
    //  The generic local search algorithm is shown below.  The various
    //  functions are overridden by the specific classes which provide
    //  the concrete details.
    //
    //  reset()         // resets curr and best
    //  while ( !isDone() ) {     // while not done
    //    beginRun()        // begin a new run
    //    do {          // do while run is not done
    //      beginStage()      // end of stage processing
    //      method = selectMethod()   // select a method
    //      switch( method ) {      // apply the method
    //        LLOYD:  curr.Lloyd(); break
    //        SWAP:   curr.Swap();  break
    //        RANDOM: curr.Random();  break
    //      }
    //      endStage()        // end of stage processing
    //    } while ( !isRunDone() )    // while run is not done
    //    tryAcceptance()     // accept if appropriate
    //    endRun()        // end of run processing
    //  }
    //  return best       // return best solution
    //
    //  Termination Parameters
    //  ----------------------
    //  Each algorithm relies on various parameter settings to determine
    //  termination conditions for runs and phases.  These settings are
    //  stored in a class KMterm (see KMterm.h).
    //
    //  Parameters (in term)
    //  --------------------
    //  maxTotStage
    //    Maximum number of stages total.
    //------------------------------------------------------------------------

    class KMeansBase {       // generic local search
    protected:
                  // fixed quantities
        int     nPts;     // number of points
        int     kCtrs;      // number of centers
        int     dim;      // dimension
        ::KMterm    term;     // termination conditions
        int     maxTotStage;    // max total stages (from term)
              // varying quantities
        int     stageNo;    // current stage number
        int     runInitStage;   // stage at which run started
        ::KMfilterCenters curr;     // current solution
        ::KMfilterCenters best;     // saved solution
    protected:          // utility functions
        virtual void printStageStats() {    // print stage information
      if (kmStatLev >= STAGE) {
                *kmOut << "\t<stage: "  << stageNo
                     << " curr: "   << curr.getAvgDist()
                     << " best: "   << best.getAvgDist()
         << " >" << endl;
      }
        }
    public:
                    // constructor
        KMeansBase(const ::KMfilterCenters &sol, const ::KMterm &t)
      : term(t), curr(sol), best(sol) {
      nPts    = sol.getNPts();
      kCtrs   = sol.getK();
      dim     = sol.getDim();
      stageNo = 0;
      maxTotStage = term.getMaxTotStage(kCtrs, nPts);
        }

        virtual ~KMeansBase() { }      // virtual destructor
        virtual ::KMfilterCenters execute();    // execute the algorithm

        int getTotalStages() const {    // return total no. of stages
          return stageNo;
        }

    protected:          // overridden by subclasses
        virtual void reset() {      // reset everythinnessg
      stageNo = 0;
      runInitStage = 0;
      curr.genRandom();     // generate random centers
      curr.getDist();       // compute initial distortion
      best = curr;
        }
        virtual bool isDone() const {   // are we done?
      return stageNo >= maxTotStage;
        }
        virtual void beginRun() {     // begin of run processing
      runInitStage = stageNo;
        }
        virtual void beginStage() { }   // start of stage processing
        virtual ::KMalg selectMethod() = 0;   // method: LLOYD or SWAP
        virtual void endStage() {     // end of stage processing
      stageNo++;
        }
        virtual bool isRunDone() {      // is run done?
      return isDone();
        }
        virtual void endRun() { }       // end of run processing
        virtual void tryAcceptance() {    // test acceptance
      if (curr.getDist() < best.getDist()) {  // is current distortion lower?
          best = curr;      // then best the current
      }
        }
    };




    //------------------------------------------------------------------------
    // KMlocalHybrid - a hybrid version combining Lloyd's and swaps.
    //
    //  Overview
    //  --------
    //  This implements a hybrid algorithm, which combines both of the
    //  previous methods along with somethinnessg like simulated annealing.
    //  Because this uses a fast annealing schedule (T = T*timeFact) it
    //  should probably be called simulated quenching.
    //
    //  The algorithm's execution is broken into the following different
    //  groups.
    //
    //  Stage:  One invocation of the Swap or Lloyd's algorithm.
    //  Run:  A run consists of a consecutive sequence of swaps
    //    followed by a consecutive sequence of Lloyd's steps.  A
    //    graphical representation of one run is presented below.
    //
    //        +--+            +---+             +--------> save ----->
    //        |  |            |   |            Y|
    //        V  |            V   |             |     N
    //  ----> Swap ------>   Lloyd's   ----> Accept? ----> restore -->
    //                                                       old
    //
    //  Decisions
    //  ---------
    //  There are three decisions which this algorithm needs to make.
    //
    //  Continue swapping or move on to Lloyd's?
    //      This is based on the simulated annealing decision choice,
    //      based on the consecutive RDL.
    //
    //  Make another pass through Lloyd's or go to acceptance?
    //      This is based on whether the relative improvement since the
    //      last stage (consecutive relative distortion loss) is above a
    //      certain fixed threshold (minConsecRDL).
    //
    //  Accept the current solution:
    //      If the current solution is better than the saved solution,
    //      then yes.  Otherwise, use the simulated annealing decision
    //      choice, based on the accumulated RDL.
    //
    //  Simulated Annealing Choice
    //  --------------------------
    //  Following the principal of simulated annealing, we somtimes
    //  chose to accept a solution, even if it is not an improvement.
    //  This choice occurs with probability:
    //
    //    exp(RDL/T),
    //
    //  where RDL is the relative distortion loss (relative to the
    //  saved solution), and T is the current temperature.  Note that
    //  if RDL > 0 (improvement) then this quantity is > 1, and so we
    //  always accept.  (Note that this is not formally correct, since
    //  in simulated annealing the probability should be based on the
    //  absolute distortion change not the relative distortion change.)
    //
    //  Annealing Schedule
    //  ------------------
    //  The temperature value T is a decreasing function of the number
    //  of the number of stages.  It starts at some initial value T0 and
    //  decreases slowly over time.  Rather than using the standard (slow)
    //  Boltzman annealing schedule, we use the following fast annealing
    //  schedule, every L stages we set
    //
    //    T = TF * T,
    //
    //  where:
    //    L (= tempRunLength) is an integer parameter set by the
    //        user.  (Presumably it depends on the number of
    //        centers and the dimension of the space.)
    //    TF (= tempReducFactor) is a real number of the form 1-x,
    //        for some small x.
    //
    //  These parameters are provided by the user and are stored in the
    //  termination condition class.
    //
    //  Initial Temperature
    //  -------------------
    //  The initial temperature T0 is a tricky thinnessg to set.  The user
    //  supplies a parameter p0 = initProbAccept, the initial
    //  probability of accepting a random swap.  However, the
    //  probability of acceting a swap depends on the given RDL value.
    //  To estimate this, for the first L runs we use p0 as the
    //  probability.  Over these runs we compute the average rdl value.
    //  Once the first L runs have completed, we set T0 so that:
    //
    //    exp(-averageRDL/T0) = initProbAccept,
    //
    //  or equivalently
    //
    //    T0 = -averageRDL/(ln initProbAccept).
    //
    //  Parameters (in term)
    //  --------------------
    //  initProbAccept
    //    Initial probability of accepting an solution that does
    //    not alter the distortion.
    //  tempRunLength
    //    The number of stages before chaning the temperature.
    //  tempReducFactor
    //    The factor by which temperature is reduced at the end of
    //    a temperature run.
    //  minConsecRDL
    //    Minimum consecutive RDL needed to keep Lloyd's algorithm
    //    going.
    //
    //  State variables
    //  ---------------
    //  temperature
    //    Temperature parameter for simulated annealing.
    //  initTempRunStage
    //    The stage number at the start of a temperature run.
    //  areSwapping
    //    True if we are still in swapping portion of a run, and
    //    false if we are in the Lloyd's portion of the run.
    //  prevDist
    //    In order to compute the consecutive RDL we need to save
    //    the distortion at the start of each trip through Lloyd's
    //    algorithm.
    //  save  Simulated annealing may move to a worse solution at
    //    times.  The "best" solution stores the global best.  The
    //    "save" solution is the solution that simulated annealing
    //    falls back to if a solution is not accepted.
    //  sumTrials
    //    Sum of rdl values during the initial trials, where we
    //    are trying to estimate the mean RDL value.
    //  trialCt
    //    A counter which starts at the number of initial trials.
    //    When it reaches 0, we compute the average RDL value and
    //    set the initial temperature value.
    //
    //  Utilities particular to Simulated Annealing
    //  -------------------------------------------
    //  simAnnealAccept()
    //    Random acceptance test for simulated annealing.
    //  initTempRun()
    //    Initialize a temperature run process by computing the
    //    initial temperature value and saving the current
    //    stage number in initTempRunStage.
    //  isTempRunDone()
    //    Check whether the current temperature run is done.
    //  endTempRun()
    //    End a temperature run by updating the temperature value
    //    and saving the current stage number in initTempRunStage.
    //
    //  Overridden methods
    //  ------------------
    //  reset()
    //    Do base class resetting.  Initialize areSwapping to
    //    true.  Save the current solution in save.  Initialize
    //    temperature runs.
    //  beginStage()
    //    Save current distortion in prevDist (for computing
    //    consecutive RDL).
    //  selectMethod()
    //    If we are swapping, use Swap, and otherwise use Lloyd's.
    //  endStage()
    //    Increment the stage number, get the distortion, and print
    //    the end-of-stage information.
    //  isRunDone()
    //    If we are swapping, then compute the simulated annealing
    //    random choice.  If it rejects, then we have ended the
    //    swapping process and are transitioning to Lloyd's
    //    algorithm.  In either case, return false.
    //
    //    If we are doing Lloyd's, return true if the consecutive
    //    RDL is less than minConsecRDL, meaning that Lloyd's
    //    algorithm has converged.
    //  endRun()
    //    If temperature run is done then do processing for the
    //    end of a temperature run.  Set areSwapping to true.
    //  tryAcceptance()
    //    At this point we have completed all the swaps and all
    //    the Lloyd's iterations.  If the distortion has improved
    //    or if the simulated annealing random choice says to,
    //    save the current solution (and update the overall best
    //    if needed).  Otherwise, restore the previously saved
    //    solution.
    //------------------------------------------------------------------------

    class KMeans : public ::KMlocal {
    private:
        double	temperature;			// temperature used in SA
        int		initTempRunStage;		// stage when temp run started
        bool	areSwapping;			// are we swapping? (or Lloyd's)
        double	prevDist;			// distortion from prev stage
        double	sumTrials;			// sum of RDL's over trials
        int		trialCt;			// trial count
        ::KMfilterCenters save;			// saved solution
    protected:					// local utilities
        double accumRDL()				// accumulated RDL
          { return (save.getDist() - curr.getDist()) /  save.getDist(); }

        double consecRDL()				// consecutive RDL
          { return (prevDist - curr.getDist()) / prevDist; }

        virtual void printStageStats() {		// print end of stage info
      if (kmStatLev >= STAGE) {
          *kmOut << "    <stage: "	<< stageNo
               << " curr: "		<< curr.getAvgDist()
               << " best: "		<< best.getAvgDist()
               << " save: "		<< save.getAvgDist()
               << " consecRDL: "	<< consecRDL()
         << " >" << endl;
      }
        }
        virtual void printRunStats() {		// print end of run info
      if (kmStatLev >= STAGE) {
          *kmOut << "    <End of Run>" << endl;
      }
        }
    protected:					// SA utilities
        int nTrials()				// number of trials
          { return kmMax(20, term.getTempRunLength()); }

        bool simAnnealAccept(double rdl) {		// random accept choice
          double prob;
          if (--trialCt >= 0) {			// still in trial phase?
            sumTrials += fabs(rdl);			// increment sum of RDLs
      if (trialCt == 0) {			// last trial?  get temp
        temperature = -sumTrials/(nTrials()*log(term.getInitProbAccept()));
        initTempRunStage = stageNo;		// start counting stages
      }
      prob = term.getInitProbAccept();	// use initial probability
          }
          else {					// use SA probability
            prob = kmMin(term.getInitProbAccept(), exp(rdl/temperature));
          }
          return prob > kmRanUnif();
        }

        void initTempRuns() {			// initialize for temp runs
          sumTrials = 0;
          trialCt = nTrials();
          initTempRunStage = ::KM_HUGE_INT;		// not counting stages
        }

        bool isTempRunDone()			// end of temperature run?
          { return stageNo - initTempRunStage >= term.getTempRunLength(); }

        void endTempRun() {				// process end of temp run
          temperature *= term.getTempReducFact();	// reduce temperature
          initTempRunStage = stageNo;
        }

    public:
                    // constructor
        KMeans(const ::KMfilterCenters &sol, const ::KMterm &t)
        : ::KMlocal(sol, t), save(sol) { }
    protected:					// overridden methods
        virtual void reset() {
          ::KMlocal::reset();				// reset base class
          save = curr;				// save initial centers
          areSwapping = true;			// start with swapping
          initTempRuns();				// initialize sim. annealing
          printStageStats();
        }
        virtual void beginStage() {			// start of stage processing
          prevDist = curr.getDist();		// save previous distortion
        }
        virtual ::KMalg selectMethod() {		// select method
          return (areSwapping ? SWAP : LLOYD );
        }
        virtual void endStage() {			// end of stage processing
          stageNo++;				// increment stage number
          curr.getDist();				// get distortion
          printStageStats();
        }
        virtual bool isRunDone() { 			// run is done
          if (areSwapping) {			// swapping?
            if (!simAnnealAccept(consecRDL())) {	// check SA acceptance
        areSwapping = false;			// transition to Lloyd's
      }
      return false;
          }
          else {					// doing Lloyd's algorithm
            return consecRDL() <= term.getMinConsecRDL();  // test for convergence
          }
        }
        virtual void endRun() { 			// end of run processing
          if (isTempRunDone()) endTempRun();	// check/process end of temp run
          areSwapping = true;			// return to swapping
          printRunStats();
        }
        virtual void tryAcceptance() { 		// test acceptance
          if (accumRDL() > 0) {			// improvement over saved?
            save = curr;				// save this one
      if (save.getDist() < best.getDist()) {	// new best?
        best = save;
      }
          }
          else if (simAnnealAccept(accumRDL())) {	// SA says save anyway
            save = best;
          }
          else {					// reject, restore old solution
            curr = save;
          }
        }
    };




  }
}

#endif
