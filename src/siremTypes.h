/*************************************************************************
 *     Copyright (C) agosto 2022 Esteban del Castillo Pérez
 *     esteban.delcastillo@urv.cat
 *
 *     This program is free software: you can redistribute it and/or modify
 *     it under the terms of the GNU General Public License as published by
 *     the Free Software Foundation, either version 3 of the License, or
 *     (at your option) any later version.
 *
 *     This program is distributed in the hope that it will be useful,
 *     but WITHOUT ANY WARRANTY; without even the implied warranty of
 *     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *     GNU General Public License for more details.
 *
 *     You should have received a copy of the GNU General Public License
 *     along with this program.  If not, see <http://www.gnu.org/licenses/>.
 **************************************************************************/
#ifndef SIREM_TYPES
#define SIREM_TYPES
//#include "types.h"

    typedef struct{
        int x, y;
    }PIXEL_XY;                  //2D pixel

  typedef struct{
        PIXEL_XY *set;
        int      size;
    }GROUP_XY;                  //pixel array


     typedef struct xGROUP{
        int 	size;
        int 	*set;
        struct xGROUP 	*group;  //group of chained integer elements
    }GROUP;

    typedef struct GROUP_F{
        int 	size;
        float 	*set;
        struct GROUP_F 	*group; //group of chained real elements
    }GROUP_F;

    typedef struct {
        int 	*set;          //items
        int 	size;          //group size.
    }LGROUP;                   //group

#endif

