//  HRSCCore: HRSC methods for Cactus
//  Copyright (C) 2014, David Radice <dradice@caltech.edu>
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


#ifndef HRSCC_METRIC_INFO_HH
#define HRSCC_METRIC_INFO_HH

namespace hrscc {

//! contains grid indices for alp, gxx, gyy and gzz
class metric_info {
    public:
        //! Cactus grid index for ADMBase::lapse
        static int idx_lapse;
        //! Cactus grid indices for the diagonal components of the metric
        static int idx_metric[3];
};

} // namespace hrscc

#endif
