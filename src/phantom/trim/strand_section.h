/*
 *  strand_section.h
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

#include "phantom/shared/segment.h"
#include "phantom/shared/strand.h"

typedef struct _strand_section {
        
        double entry_point[3], exit_point[3], *pre_point, *post_point;

        struct _strand_section *next_section, *prev_section;

        int in_roi, length_accepted, num_control_points, is_first_section;

        Segment *start_segment;
        Strand *strand;
        
} Strand_section;

int is_in_roi(double* pos, double trim_dist, int trim_to_cube);

void roi_intersect(double *intersect, double *p, double *v, double trim_dist, int trim_to_cube);

Strand_section* strand_section_alloc(Strand_section *prev_section, Strand *strand);

void strand_section_init(Strand_section *section, Strand *strand, int is_first_section);

void strand_section_free(Strand_section *section);

void set_section_entry(Strand_section *section, Segment *segment, double sphere_r,
                       int trim_to_cube);

void set_section_exit(Strand_section *section, Segment *segment, double sphere_r, int trim_to_cube);

double section_length(Strand_section *section);
