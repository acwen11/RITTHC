#!/bin/bash

cat <<EOT
//  HRSCCore: HRSC methods for Cactus
//  Copyright (C) 2011, David Radice <david.radice@aei.mpg.de>
//
//  This program is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with this program.  If not, see <http://www.gnu.org/licenses/>.


// Coefficients computed using the scripts in local/scripts/lgl.

#include <hrscc_gll_element.hh>

namespace hrscc {

EOT

for idx in `seq 0 $1`; do
    size=$[$idx+1]
    size2=$[$size*2]
    cat <<EOT
// GLLElement of degree $idx
template<>
CCTK_REAL const GLLElement<$idx>::node[$size] = {
#include "bits/nodes.$idx.inc"
};

template<>
CCTK_REAL const GLLElement<$idx>::weight[$size] = {
#include "bits/weights.$idx.inc"
};

template<>
CCTK_REAL const GLLElement<$idx>::diffop[$size][$size] = {
#include "bits/diffop.$idx.inc"
};

template<>
CCTK_REAL const GLLElement<$idx>::icoeff[$size] = {
#include "bits/icoeff.$idx.inc"
};

template<>
CCTK_REAL const GLLElement<$idx>::dltop[$size][$size] = {
#include "bits/dltop.$idx.inc"
};

template<>
CCTK_REAL const GLLElement<$idx>::idltop[$size][$size] = {
#include "bits/idltop.$idx.inc"
};

template<>
CCTK_REAL const GLLElement<$idx>::prolongation[$size2][$size] = {
#include "bits/prolongation.$idx.inc"
};

template<>
CCTK_REAL const GLLElement<$idx>::restriction[$size][$size2] = {
#include "bits/restriction.$idx.inc"
};
EOT
done

cat <<EOT
} // namespace
EOT
