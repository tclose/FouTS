/*
 *  strand_section.c
 *  Numerical Fibre Generator
 *
 *  Created by Tom Close on 25/06/08.
 *  Copyright 2008 Tom Close.
 *  Distributed under the GNU General Public Licence.
 *
 *
 *
 *  This file is part of 'Numerical Fibre Generator'.
 *
 *  'Numerical Fibre Generator' is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  'Numerical Fibre Generator' is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with 'Numerical Fibre Generator'.  If not, see <http://www.gnu.org/licenses/>
 *
 */

#include <stdlib.h>
#include <math.h>
#include <iostream>

#include "phantom/trim/strand_section.h"
#include "phantom/shared/shared.h"
#include "phantom/shared/segment.h"



Strand_section* strand_section_alloc(Strand_section *prev_section, Strand *strand) {

	Strand_section *new_section;
	
	new_section = (Strand_section*)malloc(sizeof(Strand_section));
	prev_section->next_section = new_section;
	new_section->prev_section = prev_section;
	
	strand_section_init(new_section, strand, 0);
	
	return new_section;
} 


void strand_section_init(Strand_section *section, Strand *strand, int is_first_section) {
	section->next_section = NULL;
	section->num_control_points = 0;
	section->in_roi = 0;
	section->length_accepted = 0;
	section->strand = strand;
	section->is_first_section = is_first_section;


}

void strand_section_free(Strand_section *section) {
	
	while (section->next_section != NULL) {
	
		section = section->next_section;
		free(section->prev_section);
	}
	
	free(section);

}

void set_section_entry(Strand_section *section, Segment *segment, double trim_dist, int trim_to_cube) {

  roi_intersect(section->entry_point, segment->start_point.pos, segment->disp, trim_dist, trim_to_cube);
	section->pre_point = segment->start_point.pos;
	section->in_roi = 1;
	section->start_segment = segment->next_segment;


}


int  is_in_roi(double* pos, double trim_dist, int trim_to_cube) {

  int in_roi;

  //Trim to sphere
  if (!trim_to_cube)
    in_roi = vector_norm(pos) < trim_dist;

  //Trim to cube.
  else
    in_roi = (fabs(pos[0]) < trim_dist) && (fabs(pos[1]) < trim_dist) && (fabs(pos[2]) < trim_dist);

  return in_roi;

}

void roi_intersect(double *intersect, double *p, double *v, double trim_dist, int trim_to_cube) {

  if (!trim_to_cube)
    sphere_intersect(intersect, p, v, trim_dist);

  else {

    /* NB p - point, v - vector. Finds the intersection along the line extension of 'v' about 'p' that intersects with a sphere centred on the origin with radius 'sphere_r' */

    double mu;

    if (is_in_roi(p, trim_dist, trim_to_cube)) {

      double min_mu = INFINITY;

      for (int dim_i = 0; dim_i < 3; ++dim_i) {

        int sign_v = (v[dim_i] > 0) - (v[dim_i] < 0);

        double dim_mu = (trim_dist * (double)sign_v  - p[dim_i]) / v[dim_i];

        if ((dim_mu < min_mu) && (dim_mu >= 0))
          min_mu = dim_mu;

      }

      mu = min_mu;

    } else {

      double max_mu = 0.0;

      for (int dim_i = 0; dim_i < 3; ++dim_i) {
        if (fabs(p[dim_i]) > trim_dist) {

          int sign_v = (v[dim_i] > 0) - (v[dim_i] < 0);

          double dim_mu = (- trim_dist * (double)sign_v  - p[dim_i]) / v[dim_i];

          if (dim_mu > max_mu)
            max_mu = dim_mu;

        }

      }

      mu = max_mu;

    }



    intersect[X] = p[X] + mu * v[X];
    intersect[Y] = p[Y] + mu * v[Y];
    intersect[Z] = p[Z] + mu * v[Z];

  }

}

void set_section_exit(Strand_section *section, Segment *segment, double trim_dist, int trim_to_cube) {

  roi_intersect(section->exit_point, segment->start_point.pos, segment->disp, trim_dist, trim_to_cube);
	section->post_point = segment->end_point.pos;
	section->in_roi = 0;

}

double section_length(Strand_section *section) {
	
	return dist_between_points(section->entry_point, section->exit_point);

}
