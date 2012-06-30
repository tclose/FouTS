/*
    Copyright 2009 Brain Research Institute, Melbourne, Australia

    Written by Thomas G. Close, 04/03/2009.

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

#include <list>

#include "bts/cmd.h"

#include "progressbar.h"

#include "bts/common.h"

#include "bts/fibre/strand/set.h"

#include "k_means/KMlocal.h"			// k-means algorithms




#include "bts/inline_functions.h"


using namespace BTS;

SET_VERSION_DEFAULT;
SET_AUTHOR ("Thomas G. Close");
SET_COPYRIGHT (NULL);

DESCRIPTION = {
    "Clusters strands together and returns the cluster centres as a new strand set.",
    "",
  NULL
};

ARGUMENTS = {
  Argument ("input", "The strands to be merged.").type_file (),
  Argument ("output", "The output merged strands").type_file(),
  Argument()
};


OPTIONS = {

  Option ("new_num_strands", "The number of strands ('k' in k-means) in the strands will be clustered into.")
   + Argument ("new_num_strands", "The number of strands ('k' in k-means) in the strands will be clustered into.").type_integer (1, 1000, 20),

  Option ("save_clusters", "Save the generated clusters as seperate Strand sets in the directory provided.")
   + Argument ("save_clusters", ""),

//Consult the KML documentation for the following options.

  Option ("stages", "The number of stages the k-means clustering is run for.")
   + Argument ("stages", "The number of stages the k-means clustering is run for.").type_integer (1, 1000, 20),

  Option ("max_accum_RDL","")
   + Argument ("max_accum_RDL").type_float (0.0, 1.0, 0.1),

  Option ("min_accum_RDL","")
   + Argument ("min_accum_RDL").type_float (0.0, 1.0, 0.1),

  Option ("max_run_stages","")
   + Argument ("max_run_stages").type_integer (1, 100, 3),

  Option ("init_prob_acceptance","")
   + Argument ("init_prob_acceptance").type_float (0.0, 1.0, 0.1),

  Option ("temp_run_length","")
   + Argument ("temp_run_length").type_integer (1, 100, 10),

  Option ("temp_reduc_factor","")
   + Argument ("temp_reduc_factor").type_float (0.0, 1.0, 0.95),

  Option ("num_points", "The number of points that will be generated along the strand location")
   + Argument ("points").type_integer (1, 2000, 100),

  Option ("degree", "The degree of the Strand coefficients used to describe the strands")
   + Argument ("degree").type_integer (1, 200, 3),


Option() };


void printSummary(			// print final summary
    const KMlocal&	theAlg,		// the algorithm
    const KMdata&	dataPts,	// the points
    KMfilterCenters&	ctrs);		// the centers

//void printPt(				// print a point
//    ostream&		out,		// output stream
//    const KMpoint&	p,
//    size_t dim);		// the point


EXECUTE {

  std::string input_location = argument[0];
  std::string output_location = argument[1];


  size_t new_num_strands = 20;
  bool save_clusters = false;
  std::string cluster_location = "";
  size_t num_stages = 100;
  double max_accum_RDL = 0.1;
  double min_accum_RDL = 0.1;
  int max_run_stages = 3;
  double init_prob_acceptance = 0.5;
  int temp_run_length = 10;
  double temp_reduc_factor = 0.95;
  size_t degree = 0;
  size_t num_points = 0;

  Options opt = get_options("new_num_strands");
  if (opt.size())
    new_num_strands = opt[0][0];

  opt = get_options("cluster_location");
  if (opt.size()) {
    save_clusters = true;
    cluster_location = opt[0][0].c_str();
  }

  opt = get_options("num_stages");
  if (opt.size())
    num_stages = opt[0][0];

  opt = get_options("max_accum_RDL");
  if (opt.size())
    max_accum_RDL = opt[0][0];

  opt = get_options("min_accum_RDL");
  if (opt.size())
    min_accum_RDL = opt[0][0];

  opt = get_options("max_run_stages");
  if (opt.size())
    max_run_stages = opt[0][0];

  opt = get_options("init_prob_acceptance");
  if (opt.size())
    init_prob_acceptance = opt[0][0];

  opt = get_options("temp_run_length");
  if (opt.size())
    temp_run_length = opt[0][0];

  opt = get_options("temp_reduc_factor");
  if (opt.size())
    temp_reduc_factor = opt[0][0];

 opt = get_options("num_points");
  if (opt.size())
    num_points = opt[0][0];

  opt = get_options("degree");
  if (opt.size())
    degree = opt[0][0];


  Fibre::Strand::Set input_strands(input_location, degree);

  if (new_num_strands >= input_strands.size())
    throw Exception("New number of strands (" + str(new_num_strands) + ") must be less than the original number of strands (" + str(input_strands.size()) + ")");


  KMterm	term(num_stages, 0, 0, 0,		// run for 100 stages
      min_accum_RDL,			// min consec RDL
      max_accum_RDL,			// min accum RDL
      max_run_stages,			// max run stages
      init_prob_acceptance,			// init. prob. of acceptance
      temp_run_length,			// temp. run length
      temp_reduc_factor);			// temp. reduction factor



//-------------------------//
//  k-means cluster strands
//-------------------------//

  degree = input_strands[0].degree();

  int	ptDim	= degree * 3;		// dimension

  KMdata dataPts(ptDim, input_strands.size());		// allocate data storage

  for (size_t strand_i = 0; strand_i < input_strands.size(); strand_i++) {
    for (size_t degree_i = 0; degree_i < input_strands[strand_i].degree(); degree_i++) {
      dataPts[strand_i][degree_i * 3 + X] = input_strands[strand_i][degree_i][X];
      dataPts[strand_i][degree_i * 3 + Y] = input_strands[strand_i][degree_i][Y];
      dataPts[strand_i][degree_i * 3 + Z] = input_strands[strand_i][degree_i][Z];
    }
  }



  dataPts.buildKcTree();			// build filtering structure

  KMfilterCenters ctrs(new_num_strands, dataPts);		// allocate centers
//
//  KMlocalLloyds kmLloyds(ctrs, term);		// repeated Lloyd's
//  ctrs = kmLloyds.execute();			// execute
//
//  KMlocalSwap kmSwap(ctrs, term);		// Swap heuristic
//  ctrs = kmSwap.execute();
//
//  KMlocalEZ_Hybrid kmEZ_Hybrid(ctrs, term);	// EZ-Hybrid heuristic
//  ctrs = kmEZ_Hybrid.execute();

  KMlocalHybrid kmHybrid(ctrs, term);		// Hybrid heuristic
  ctrs = kmHybrid.execute();

//--------------------------//
//  Get cluster information
//--------------------------//


  KMctrIdxArray closeCtr = new KMctrIdx[dataPts.getNPts()];
  double* sqDist = new double[dataPts.getNPts()];
  ctrs.getAssignments(closeCtr, sqDist);




//------------------------------//
//  Print results of clustering
//------------------------------//

  for (int ctr_i = 0; ctr_i < ctrs.getK(); ctr_i++) {

    std::cout << "Centre: " << ctr_i << " [";

    for (int pt_i = 0; pt_i < dataPts.getNPts(); pt_i++) {
      if (closeCtr[pt_i] == ctr_i)
        std::cout << pt_i << ",";
    }

    std::cout << "] - [";

    for (int pt_i = 0; pt_i < dataPts.getNPts(); pt_i++) {
      if (closeCtr[pt_i] == ctr_i)
        std::cout << setprecision(2) << sqDist[pt_i] << ",";
    }

    std::cout << "]";

    std::cout << std::endl;
  }


//-----------------------------------//
//  Save centre points as new strands
//-----------------------------------//

  Fibre::Strand::Set output_strands(input_strands.get_extend_props());

  for (size_t ctr_i = 0; ctr_i < new_num_strands; ctr_i++) {

    Fibre::Strand strand(degree);

    for (size_t degree_i = 0; degree_i < input_strands[0].degree(); degree_i++) {
      for (size_t dim_i = 0; dim_i < 3; dim_i++)
        strand[degree_i][dim_i] = ctrs[ctr_i][degree_i * 3 + dim_i];
    }

    output_strands.push_back(strand);

  }

  output_strands.save(output_location, num_points);



//----------------------------------------------//
//  Save clusters in separate files if required
//----------------------------------------------//


  if (save_clusters) {

    size_t num_cluster_dec_places = num_dec_places(new_num_strands);

    std::vector<Fibre::Strand::Set> clusters;

    for (size_t ctr_i = 0; ctr_i < new_num_strands; ctr_i++)
      clusters.push_back(Fibre::Strand::Set());


    for (size_t strand_i = 0; strand_i < input_strands.size(); strand_i++)
      clusters[closeCtr[strand_i]].push_back(input_strands[strand_i]);


    for (size_t ctr_i = 0; ctr_i < new_num_strands; ctr_i++)
      clusters[ctr_i].save(File::join(cluster_location, "cluster_" + str(ctr_i, num_cluster_dec_places) + "." + Fibre::Strand::FILE_EXTENSION));


  }

  delete [] closeCtr;
  delete [] sqDist;


}







void printPt(ostream& out, const KMpoint& p, int dim)
{
    out << "(" << p[0];
    for (int i = 1; i < dim; i++) {
	out << ", " << p[i];
    }
    out << ")\n";
}

//------------------------------------------------------------------------
//  Print summary of execution
//------------------------------------------------------------------------
void printSummary(
    const KMlocal&		theAlg,		// the algorithm
    const KMdata&		dataPts,	// the points
    KMfilterCenters&		ctrs)		// the centers
{
    cout << "Number of stages: " << theAlg.getTotalStages() << "\n";
    cout << "Average distortion: " <<
	         ctrs.getDist(false)/double(ctrs.getNPts()) << "\n";
    					// print final center points
    cout << "(Final Center Points:\n";
    ctrs.print();
    cout << ")\n";
    					// get/print final cluster assignments
    KMctrIdxArray closeCtr = new KMctrIdx[dataPts.getNPts()];
    double* sqDist = new double[dataPts.getNPts()];
    ctrs.getAssignments(closeCtr, sqDist);

    *kmOut	<< "(Cluster assignments:\n"
		<< "    Point  Center  Squared Dist\n"
		<< "    -----  ------  ------------\n";
    for (int i = 0; i < dataPts.getNPts(); i++) {
	*kmOut	<< "   " << setw(5) << i
		<< "   " << setw(5) << closeCtr[i]
		<< "   " << setw(10) << sqDist[i]
		<< "\n";
    }
    *kmOut << ")\n";
    delete [] closeCtr;
    delete [] sqDist;
}
