/*
 * A mex wrapper for the mrtrix load tracks function.
 
 Can be compiled with the following line:
 
 mex load_mrtrix_strands.cpp -I/data/home/tclose/Code/Tractography/mrtrix_bundle/src/mrtrix/lib -I/data/home/tclose/Code/Tractography/mrtrix_bundle/include -I/data/home/tclose/Code/Tractography/mrtrix_bundle/src/mrtrix/src -I/data/home/tclose/Code/Tractography/tractography/src -DMACOSX -DBYTE_ORDER_IS_BIG_ENDIAN -DDEBUG_EXTENSIONS
 
 */
 
#include <iostream>
#include <math.h>
#include "mex.h"
#include "app.h"
#include "math/matrix.h"
#include "math/least_squares.h"
#include "point.h"
#include "file/path.h"
#include "progressbar.h"
#include "shared.h"



extern void _main();


static int count_mrtrix_strands(char *filename) {
  MR::DWI::Tractography::Properties properties;
  std::vector<MR::Point> tck;
  
  MR::DWI::Tractography::Reader<double> reader;
  
  reader.open(filename, properties);
    
  int count = 0;
  
  while (reader.next(tck)) {
   
    if (tck.size()) {
      count++;

    }
    
  }
        
  reader.close ();
  
  return count;
}


void mexFunction(
		 int          nlhs,
		 mxArray      *lhs[],
		 int          nrhs,
		 const mxArray *prhs[]
		 )
{
  double      *vin1, *vin2;

  /* Check for proper number of arguments */

  if (nrhs != 1) {
    mexErrMsgTxt("MEXCPP requires one input argument (filename).");
  } 

  char filename[500];

  if (mxGetString(prhs[0], filename, 499)) {
    mexErrMsgTxt("Input argument was not valid string.");
  }

  int strand_count = count_mrtrix_strands(filename);
  
  lhs[0] = mxCreateCellMatrix(strand_count, 1);
  
  MR::DWI::Tractography::Properties properties;
  std::vector<MR::Point> tck;
  
  MR::DWI::Tractography::Reader<double> reader;
  
  reader.open(filename, properties);
    
  int count = 0;
  
  while (reader.next(tck)) {
   
    if (tck.size()) {
      
      mxArray *tck_mxArray = mxCreateDoubleMatrix(3,(tck.size()), mxREAL);
      
      double *tck_matrix = mxGetPr(tck_mxArray);
      
      
      for (size_t point_i = 0; point_i < tck.size(); point_i++) {
        tck_matrix[point_i*3 + X] = tck[point_i][X];
        tck_matrix[point_i*3 + Y] = tck[point_i][Y];
        tck_matrix[point_i*3 + Z] = tck[point_i][Z];              
      }

      mxSetCell(lhs[0], count, tck_mxArray);

      count++;

    }
    
  }
        
  reader.close ();
  

  
  
  
  
  return;
}
