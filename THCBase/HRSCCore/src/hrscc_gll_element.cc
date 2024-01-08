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

// GLLElement of degree 0
template<>
CCTK_REAL const GLLElement<0>::node[1] = {
#include "bits/nodes.0.inc"
};

template<>
CCTK_REAL const GLLElement<0>::weight[1] = {
#include "bits/weights.0.inc"
};

template<>
CCTK_REAL const GLLElement<0>::diffop[1][1] = {
#include "bits/diffop.0.inc"
};

template<>
CCTK_REAL const GLLElement<0>::icoeff[1] = {
#include "bits/icoeff.0.inc"
};

template<>
CCTK_REAL const GLLElement<0>::dltop[1][1] = {
#include "bits/dltop.0.inc"
};

template<>
CCTK_REAL const GLLElement<0>::idltop[1][1] = {
#include "bits/idltop.0.inc"
};

template<>
CCTK_REAL const GLLElement<0>::prolongation[2][1] = {
#include "bits/prolongation.0.inc"
};

template<>
CCTK_REAL const GLLElement<0>::restriction[1][2] = {
#include "bits/restriction.0.inc"
};
// GLLElement of degree 1
template<>
CCTK_REAL const GLLElement<1>::node[2] = {
#include "bits/nodes.1.inc"
};

template<>
CCTK_REAL const GLLElement<1>::weight[2] = {
#include "bits/weights.1.inc"
};

template<>
CCTK_REAL const GLLElement<1>::diffop[2][2] = {
#include "bits/diffop.1.inc"
};

template<>
CCTK_REAL const GLLElement<1>::icoeff[2] = {
#include "bits/icoeff.1.inc"
};

template<>
CCTK_REAL const GLLElement<1>::dltop[2][2] = {
#include "bits/dltop.1.inc"
};

template<>
CCTK_REAL const GLLElement<1>::idltop[2][2] = {
#include "bits/idltop.1.inc"
};

template<>
CCTK_REAL const GLLElement<1>::prolongation[4][2] = {
#include "bits/prolongation.1.inc"
};

template<>
CCTK_REAL const GLLElement<1>::restriction[2][4] = {
#include "bits/restriction.1.inc"
};
// GLLElement of degree 2
template<>
CCTK_REAL const GLLElement<2>::node[3] = {
#include "bits/nodes.2.inc"
};

template<>
CCTK_REAL const GLLElement<2>::weight[3] = {
#include "bits/weights.2.inc"
};

template<>
CCTK_REAL const GLLElement<2>::diffop[3][3] = {
#include "bits/diffop.2.inc"
};

template<>
CCTK_REAL const GLLElement<2>::icoeff[3] = {
#include "bits/icoeff.2.inc"
};

template<>
CCTK_REAL const GLLElement<2>::dltop[3][3] = {
#include "bits/dltop.2.inc"
};

template<>
CCTK_REAL const GLLElement<2>::idltop[3][3] = {
#include "bits/idltop.2.inc"
};

template<>
CCTK_REAL const GLLElement<2>::prolongation[6][3] = {
#include "bits/prolongation.2.inc"
};

template<>
CCTK_REAL const GLLElement<2>::restriction[3][6] = {
#include "bits/restriction.2.inc"
};
// GLLElement of degree 3
template<>
CCTK_REAL const GLLElement<3>::node[4] = {
#include "bits/nodes.3.inc"
};

template<>
CCTK_REAL const GLLElement<3>::weight[4] = {
#include "bits/weights.3.inc"
};

template<>
CCTK_REAL const GLLElement<3>::diffop[4][4] = {
#include "bits/diffop.3.inc"
};

template<>
CCTK_REAL const GLLElement<3>::icoeff[4] = {
#include "bits/icoeff.3.inc"
};

template<>
CCTK_REAL const GLLElement<3>::dltop[4][4] = {
#include "bits/dltop.3.inc"
};

template<>
CCTK_REAL const GLLElement<3>::idltop[4][4] = {
#include "bits/idltop.3.inc"
};

template<>
CCTK_REAL const GLLElement<3>::prolongation[8][4] = {
#include "bits/prolongation.3.inc"
};

template<>
CCTK_REAL const GLLElement<3>::restriction[4][8] = {
#include "bits/restriction.3.inc"
};
// GLLElement of degree 4
template<>
CCTK_REAL const GLLElement<4>::node[5] = {
#include "bits/nodes.4.inc"
};

template<>
CCTK_REAL const GLLElement<4>::weight[5] = {
#include "bits/weights.4.inc"
};

template<>
CCTK_REAL const GLLElement<4>::diffop[5][5] = {
#include "bits/diffop.4.inc"
};

template<>
CCTK_REAL const GLLElement<4>::icoeff[5] = {
#include "bits/icoeff.4.inc"
};

template<>
CCTK_REAL const GLLElement<4>::dltop[5][5] = {
#include "bits/dltop.4.inc"
};

template<>
CCTK_REAL const GLLElement<4>::idltop[5][5] = {
#include "bits/idltop.4.inc"
};

template<>
CCTK_REAL const GLLElement<4>::prolongation[10][5] = {
#include "bits/prolongation.4.inc"
};

template<>
CCTK_REAL const GLLElement<4>::restriction[5][10] = {
#include "bits/restriction.4.inc"
};
// GLLElement of degree 5
template<>
CCTK_REAL const GLLElement<5>::node[6] = {
#include "bits/nodes.5.inc"
};

template<>
CCTK_REAL const GLLElement<5>::weight[6] = {
#include "bits/weights.5.inc"
};

template<>
CCTK_REAL const GLLElement<5>::diffop[6][6] = {
#include "bits/diffop.5.inc"
};

template<>
CCTK_REAL const GLLElement<5>::icoeff[6] = {
#include "bits/icoeff.5.inc"
};

template<>
CCTK_REAL const GLLElement<5>::dltop[6][6] = {
#include "bits/dltop.5.inc"
};

template<>
CCTK_REAL const GLLElement<5>::idltop[6][6] = {
#include "bits/idltop.5.inc"
};

template<>
CCTK_REAL const GLLElement<5>::prolongation[12][6] = {
#include "bits/prolongation.5.inc"
};

template<>
CCTK_REAL const GLLElement<5>::restriction[6][12] = {
#include "bits/restriction.5.inc"
};
// GLLElement of degree 6
template<>
CCTK_REAL const GLLElement<6>::node[7] = {
#include "bits/nodes.6.inc"
};

template<>
CCTK_REAL const GLLElement<6>::weight[7] = {
#include "bits/weights.6.inc"
};

template<>
CCTK_REAL const GLLElement<6>::diffop[7][7] = {
#include "bits/diffop.6.inc"
};

template<>
CCTK_REAL const GLLElement<6>::icoeff[7] = {
#include "bits/icoeff.6.inc"
};

template<>
CCTK_REAL const GLLElement<6>::dltop[7][7] = {
#include "bits/dltop.6.inc"
};

template<>
CCTK_REAL const GLLElement<6>::idltop[7][7] = {
#include "bits/idltop.6.inc"
};

template<>
CCTK_REAL const GLLElement<6>::prolongation[14][7] = {
#include "bits/prolongation.6.inc"
};

template<>
CCTK_REAL const GLLElement<6>::restriction[7][14] = {
#include "bits/restriction.6.inc"
};
// GLLElement of degree 7
template<>
CCTK_REAL const GLLElement<7>::node[8] = {
#include "bits/nodes.7.inc"
};

template<>
CCTK_REAL const GLLElement<7>::weight[8] = {
#include "bits/weights.7.inc"
};

template<>
CCTK_REAL const GLLElement<7>::diffop[8][8] = {
#include "bits/diffop.7.inc"
};

template<>
CCTK_REAL const GLLElement<7>::icoeff[8] = {
#include "bits/icoeff.7.inc"
};

template<>
CCTK_REAL const GLLElement<7>::dltop[8][8] = {
#include "bits/dltop.7.inc"
};

template<>
CCTK_REAL const GLLElement<7>::idltop[8][8] = {
#include "bits/idltop.7.inc"
};

template<>
CCTK_REAL const GLLElement<7>::prolongation[16][8] = {
#include "bits/prolongation.7.inc"
};

template<>
CCTK_REAL const GLLElement<7>::restriction[8][16] = {
#include "bits/restriction.7.inc"
};
// GLLElement of degree 8
template<>
CCTK_REAL const GLLElement<8>::node[9] = {
#include "bits/nodes.8.inc"
};

template<>
CCTK_REAL const GLLElement<8>::weight[9] = {
#include "bits/weights.8.inc"
};

template<>
CCTK_REAL const GLLElement<8>::diffop[9][9] = {
#include "bits/diffop.8.inc"
};

template<>
CCTK_REAL const GLLElement<8>::icoeff[9] = {
#include "bits/icoeff.8.inc"
};

template<>
CCTK_REAL const GLLElement<8>::dltop[9][9] = {
#include "bits/dltop.8.inc"
};

template<>
CCTK_REAL const GLLElement<8>::idltop[9][9] = {
#include "bits/idltop.8.inc"
};

template<>
CCTK_REAL const GLLElement<8>::prolongation[18][9] = {
#include "bits/prolongation.8.inc"
};

template<>
CCTK_REAL const GLLElement<8>::restriction[9][18] = {
#include "bits/restriction.8.inc"
};
// GLLElement of degree 9
template<>
CCTK_REAL const GLLElement<9>::node[10] = {
#include "bits/nodes.9.inc"
};

template<>
CCTK_REAL const GLLElement<9>::weight[10] = {
#include "bits/weights.9.inc"
};

template<>
CCTK_REAL const GLLElement<9>::diffop[10][10] = {
#include "bits/diffop.9.inc"
};

template<>
CCTK_REAL const GLLElement<9>::icoeff[10] = {
#include "bits/icoeff.9.inc"
};

template<>
CCTK_REAL const GLLElement<9>::dltop[10][10] = {
#include "bits/dltop.9.inc"
};

template<>
CCTK_REAL const GLLElement<9>::idltop[10][10] = {
#include "bits/idltop.9.inc"
};

template<>
CCTK_REAL const GLLElement<9>::prolongation[20][10] = {
#include "bits/prolongation.9.inc"
};

template<>
CCTK_REAL const GLLElement<9>::restriction[10][20] = {
#include "bits/restriction.9.inc"
};
// GLLElement of degree 10
template<>
CCTK_REAL const GLLElement<10>::node[11] = {
#include "bits/nodes.10.inc"
};

template<>
CCTK_REAL const GLLElement<10>::weight[11] = {
#include "bits/weights.10.inc"
};

template<>
CCTK_REAL const GLLElement<10>::diffop[11][11] = {
#include "bits/diffop.10.inc"
};

template<>
CCTK_REAL const GLLElement<10>::icoeff[11] = {
#include "bits/icoeff.10.inc"
};

template<>
CCTK_REAL const GLLElement<10>::dltop[11][11] = {
#include "bits/dltop.10.inc"
};

template<>
CCTK_REAL const GLLElement<10>::idltop[11][11] = {
#include "bits/idltop.10.inc"
};

template<>
CCTK_REAL const GLLElement<10>::prolongation[22][11] = {
#include "bits/prolongation.10.inc"
};

template<>
CCTK_REAL const GLLElement<10>::restriction[11][22] = {
#include "bits/restriction.10.inc"
};
// GLLElement of degree 11
template<>
CCTK_REAL const GLLElement<11>::node[12] = {
#include "bits/nodes.11.inc"
};

template<>
CCTK_REAL const GLLElement<11>::weight[12] = {
#include "bits/weights.11.inc"
};

template<>
CCTK_REAL const GLLElement<11>::diffop[12][12] = {
#include "bits/diffop.11.inc"
};

template<>
CCTK_REAL const GLLElement<11>::icoeff[12] = {
#include "bits/icoeff.11.inc"
};

template<>
CCTK_REAL const GLLElement<11>::dltop[12][12] = {
#include "bits/dltop.11.inc"
};

template<>
CCTK_REAL const GLLElement<11>::idltop[12][12] = {
#include "bits/idltop.11.inc"
};

template<>
CCTK_REAL const GLLElement<11>::prolongation[24][12] = {
#include "bits/prolongation.11.inc"
};

template<>
CCTK_REAL const GLLElement<11>::restriction[12][24] = {
#include "bits/restriction.11.inc"
};
// GLLElement of degree 12
template<>
CCTK_REAL const GLLElement<12>::node[13] = {
#include "bits/nodes.12.inc"
};

template<>
CCTK_REAL const GLLElement<12>::weight[13] = {
#include "bits/weights.12.inc"
};

template<>
CCTK_REAL const GLLElement<12>::diffop[13][13] = {
#include "bits/diffop.12.inc"
};

template<>
CCTK_REAL const GLLElement<12>::icoeff[13] = {
#include "bits/icoeff.12.inc"
};

template<>
CCTK_REAL const GLLElement<12>::dltop[13][13] = {
#include "bits/dltop.12.inc"
};

template<>
CCTK_REAL const GLLElement<12>::idltop[13][13] = {
#include "bits/idltop.12.inc"
};

template<>
CCTK_REAL const GLLElement<12>::prolongation[26][13] = {
#include "bits/prolongation.12.inc"
};

template<>
CCTK_REAL const GLLElement<12>::restriction[13][26] = {
#include "bits/restriction.12.inc"
};
// GLLElement of degree 13
template<>
CCTK_REAL const GLLElement<13>::node[14] = {
#include "bits/nodes.13.inc"
};

template<>
CCTK_REAL const GLLElement<13>::weight[14] = {
#include "bits/weights.13.inc"
};

template<>
CCTK_REAL const GLLElement<13>::diffop[14][14] = {
#include "bits/diffop.13.inc"
};

template<>
CCTK_REAL const GLLElement<13>::icoeff[14] = {
#include "bits/icoeff.13.inc"
};

template<>
CCTK_REAL const GLLElement<13>::dltop[14][14] = {
#include "bits/dltop.13.inc"
};

template<>
CCTK_REAL const GLLElement<13>::idltop[14][14] = {
#include "bits/idltop.13.inc"
};

template<>
CCTK_REAL const GLLElement<13>::prolongation[28][14] = {
#include "bits/prolongation.13.inc"
};

template<>
CCTK_REAL const GLLElement<13>::restriction[14][28] = {
#include "bits/restriction.13.inc"
};
// GLLElement of degree 14
template<>
CCTK_REAL const GLLElement<14>::node[15] = {
#include "bits/nodes.14.inc"
};

template<>
CCTK_REAL const GLLElement<14>::weight[15] = {
#include "bits/weights.14.inc"
};

template<>
CCTK_REAL const GLLElement<14>::diffop[15][15] = {
#include "bits/diffop.14.inc"
};

template<>
CCTK_REAL const GLLElement<14>::icoeff[15] = {
#include "bits/icoeff.14.inc"
};

template<>
CCTK_REAL const GLLElement<14>::dltop[15][15] = {
#include "bits/dltop.14.inc"
};

template<>
CCTK_REAL const GLLElement<14>::idltop[15][15] = {
#include "bits/idltop.14.inc"
};

template<>
CCTK_REAL const GLLElement<14>::prolongation[30][15] = {
#include "bits/prolongation.14.inc"
};

template<>
CCTK_REAL const GLLElement<14>::restriction[15][30] = {
#include "bits/restriction.14.inc"
};
// GLLElement of degree 15
template<>
CCTK_REAL const GLLElement<15>::node[16] = {
#include "bits/nodes.15.inc"
};

template<>
CCTK_REAL const GLLElement<15>::weight[16] = {
#include "bits/weights.15.inc"
};

template<>
CCTK_REAL const GLLElement<15>::diffop[16][16] = {
#include "bits/diffop.15.inc"
};

template<>
CCTK_REAL const GLLElement<15>::icoeff[16] = {
#include "bits/icoeff.15.inc"
};

template<>
CCTK_REAL const GLLElement<15>::dltop[16][16] = {
#include "bits/dltop.15.inc"
};

template<>
CCTK_REAL const GLLElement<15>::idltop[16][16] = {
#include "bits/idltop.15.inc"
};

template<>
CCTK_REAL const GLLElement<15>::prolongation[32][16] = {
#include "bits/prolongation.15.inc"
};

template<>
CCTK_REAL const GLLElement<15>::restriction[16][32] = {
#include "bits/restriction.15.inc"
};
// GLLElement of degree 16
template<>
CCTK_REAL const GLLElement<16>::node[17] = {
#include "bits/nodes.16.inc"
};

template<>
CCTK_REAL const GLLElement<16>::weight[17] = {
#include "bits/weights.16.inc"
};

template<>
CCTK_REAL const GLLElement<16>::diffop[17][17] = {
#include "bits/diffop.16.inc"
};

template<>
CCTK_REAL const GLLElement<16>::icoeff[17] = {
#include "bits/icoeff.16.inc"
};

template<>
CCTK_REAL const GLLElement<16>::dltop[17][17] = {
#include "bits/dltop.16.inc"
};

template<>
CCTK_REAL const GLLElement<16>::idltop[17][17] = {
#include "bits/idltop.16.inc"
};

template<>
CCTK_REAL const GLLElement<16>::prolongation[34][17] = {
#include "bits/prolongation.16.inc"
};

template<>
CCTK_REAL const GLLElement<16>::restriction[17][34] = {
#include "bits/restriction.16.inc"
};
// GLLElement of degree 17
template<>
CCTK_REAL const GLLElement<17>::node[18] = {
#include "bits/nodes.17.inc"
};

template<>
CCTK_REAL const GLLElement<17>::weight[18] = {
#include "bits/weights.17.inc"
};

template<>
CCTK_REAL const GLLElement<17>::diffop[18][18] = {
#include "bits/diffop.17.inc"
};

template<>
CCTK_REAL const GLLElement<17>::icoeff[18] = {
#include "bits/icoeff.17.inc"
};

template<>
CCTK_REAL const GLLElement<17>::dltop[18][18] = {
#include "bits/dltop.17.inc"
};

template<>
CCTK_REAL const GLLElement<17>::idltop[18][18] = {
#include "bits/idltop.17.inc"
};

template<>
CCTK_REAL const GLLElement<17>::prolongation[36][18] = {
#include "bits/prolongation.17.inc"
};

template<>
CCTK_REAL const GLLElement<17>::restriction[18][36] = {
#include "bits/restriction.17.inc"
};
// GLLElement of degree 18
template<>
CCTK_REAL const GLLElement<18>::node[19] = {
#include "bits/nodes.18.inc"
};

template<>
CCTK_REAL const GLLElement<18>::weight[19] = {
#include "bits/weights.18.inc"
};

template<>
CCTK_REAL const GLLElement<18>::diffop[19][19] = {
#include "bits/diffop.18.inc"
};

template<>
CCTK_REAL const GLLElement<18>::icoeff[19] = {
#include "bits/icoeff.18.inc"
};

template<>
CCTK_REAL const GLLElement<18>::dltop[19][19] = {
#include "bits/dltop.18.inc"
};

template<>
CCTK_REAL const GLLElement<18>::idltop[19][19] = {
#include "bits/idltop.18.inc"
};

template<>
CCTK_REAL const GLLElement<18>::prolongation[38][19] = {
#include "bits/prolongation.18.inc"
};

template<>
CCTK_REAL const GLLElement<18>::restriction[19][38] = {
#include "bits/restriction.18.inc"
};
// GLLElement of degree 19
template<>
CCTK_REAL const GLLElement<19>::node[20] = {
#include "bits/nodes.19.inc"
};

template<>
CCTK_REAL const GLLElement<19>::weight[20] = {
#include "bits/weights.19.inc"
};

template<>
CCTK_REAL const GLLElement<19>::diffop[20][20] = {
#include "bits/diffop.19.inc"
};

template<>
CCTK_REAL const GLLElement<19>::icoeff[20] = {
#include "bits/icoeff.19.inc"
};

template<>
CCTK_REAL const GLLElement<19>::dltop[20][20] = {
#include "bits/dltop.19.inc"
};

template<>
CCTK_REAL const GLLElement<19>::idltop[20][20] = {
#include "bits/idltop.19.inc"
};

template<>
CCTK_REAL const GLLElement<19>::prolongation[40][20] = {
#include "bits/prolongation.19.inc"
};

template<>
CCTK_REAL const GLLElement<19>::restriction[20][40] = {
#include "bits/restriction.19.inc"
};
// GLLElement of degree 20
template<>
CCTK_REAL const GLLElement<20>::node[21] = {
#include "bits/nodes.20.inc"
};

template<>
CCTK_REAL const GLLElement<20>::weight[21] = {
#include "bits/weights.20.inc"
};

template<>
CCTK_REAL const GLLElement<20>::diffop[21][21] = {
#include "bits/diffop.20.inc"
};

template<>
CCTK_REAL const GLLElement<20>::icoeff[21] = {
#include "bits/icoeff.20.inc"
};

template<>
CCTK_REAL const GLLElement<20>::dltop[21][21] = {
#include "bits/dltop.20.inc"
};

template<>
CCTK_REAL const GLLElement<20>::idltop[21][21] = {
#include "bits/idltop.20.inc"
};

template<>
CCTK_REAL const GLLElement<20>::prolongation[42][21] = {
#include "bits/prolongation.20.inc"
};

template<>
CCTK_REAL const GLLElement<20>::restriction[21][42] = {
#include "bits/restriction.20.inc"
};
// GLLElement of degree 21
template<>
CCTK_REAL const GLLElement<21>::node[22] = {
#include "bits/nodes.21.inc"
};

template<>
CCTK_REAL const GLLElement<21>::weight[22] = {
#include "bits/weights.21.inc"
};

template<>
CCTK_REAL const GLLElement<21>::diffop[22][22] = {
#include "bits/diffop.21.inc"
};

template<>
CCTK_REAL const GLLElement<21>::icoeff[22] = {
#include "bits/icoeff.21.inc"
};

template<>
CCTK_REAL const GLLElement<21>::dltop[22][22] = {
#include "bits/dltop.21.inc"
};

template<>
CCTK_REAL const GLLElement<21>::idltop[22][22] = {
#include "bits/idltop.21.inc"
};

template<>
CCTK_REAL const GLLElement<21>::prolongation[44][22] = {
#include "bits/prolongation.21.inc"
};

template<>
CCTK_REAL const GLLElement<21>::restriction[22][44] = {
#include "bits/restriction.21.inc"
};
// GLLElement of degree 22
template<>
CCTK_REAL const GLLElement<22>::node[23] = {
#include "bits/nodes.22.inc"
};

template<>
CCTK_REAL const GLLElement<22>::weight[23] = {
#include "bits/weights.22.inc"
};

template<>
CCTK_REAL const GLLElement<22>::diffop[23][23] = {
#include "bits/diffop.22.inc"
};

template<>
CCTK_REAL const GLLElement<22>::icoeff[23] = {
#include "bits/icoeff.22.inc"
};

template<>
CCTK_REAL const GLLElement<22>::dltop[23][23] = {
#include "bits/dltop.22.inc"
};

template<>
CCTK_REAL const GLLElement<22>::idltop[23][23] = {
#include "bits/idltop.22.inc"
};

template<>
CCTK_REAL const GLLElement<22>::prolongation[46][23] = {
#include "bits/prolongation.22.inc"
};

template<>
CCTK_REAL const GLLElement<22>::restriction[23][46] = {
#include "bits/restriction.22.inc"
};
// GLLElement of degree 23
template<>
CCTK_REAL const GLLElement<23>::node[24] = {
#include "bits/nodes.23.inc"
};

template<>
CCTK_REAL const GLLElement<23>::weight[24] = {
#include "bits/weights.23.inc"
};

template<>
CCTK_REAL const GLLElement<23>::diffop[24][24] = {
#include "bits/diffop.23.inc"
};

template<>
CCTK_REAL const GLLElement<23>::icoeff[24] = {
#include "bits/icoeff.23.inc"
};

template<>
CCTK_REAL const GLLElement<23>::dltop[24][24] = {
#include "bits/dltop.23.inc"
};

template<>
CCTK_REAL const GLLElement<23>::idltop[24][24] = {
#include "bits/idltop.23.inc"
};

template<>
CCTK_REAL const GLLElement<23>::prolongation[48][24] = {
#include "bits/prolongation.23.inc"
};

template<>
CCTK_REAL const GLLElement<23>::restriction[24][48] = {
#include "bits/restriction.23.inc"
};
// GLLElement of degree 24
template<>
CCTK_REAL const GLLElement<24>::node[25] = {
#include "bits/nodes.24.inc"
};

template<>
CCTK_REAL const GLLElement<24>::weight[25] = {
#include "bits/weights.24.inc"
};

template<>
CCTK_REAL const GLLElement<24>::diffop[25][25] = {
#include "bits/diffop.24.inc"
};

template<>
CCTK_REAL const GLLElement<24>::icoeff[25] = {
#include "bits/icoeff.24.inc"
};

template<>
CCTK_REAL const GLLElement<24>::dltop[25][25] = {
#include "bits/dltop.24.inc"
};

template<>
CCTK_REAL const GLLElement<24>::idltop[25][25] = {
#include "bits/idltop.24.inc"
};

template<>
CCTK_REAL const GLLElement<24>::prolongation[50][25] = {
#include "bits/prolongation.24.inc"
};

template<>
CCTK_REAL const GLLElement<24>::restriction[25][50] = {
#include "bits/restriction.24.inc"
};
// GLLElement of degree 25
template<>
CCTK_REAL const GLLElement<25>::node[26] = {
#include "bits/nodes.25.inc"
};

template<>
CCTK_REAL const GLLElement<25>::weight[26] = {
#include "bits/weights.25.inc"
};

template<>
CCTK_REAL const GLLElement<25>::diffop[26][26] = {
#include "bits/diffop.25.inc"
};

template<>
CCTK_REAL const GLLElement<25>::icoeff[26] = {
#include "bits/icoeff.25.inc"
};

template<>
CCTK_REAL const GLLElement<25>::dltop[26][26] = {
#include "bits/dltop.25.inc"
};

template<>
CCTK_REAL const GLLElement<25>::idltop[26][26] = {
#include "bits/idltop.25.inc"
};

template<>
CCTK_REAL const GLLElement<25>::prolongation[52][26] = {
#include "bits/prolongation.25.inc"
};

template<>
CCTK_REAL const GLLElement<25>::restriction[26][52] = {
#include "bits/restriction.25.inc"
};
} // namespace
