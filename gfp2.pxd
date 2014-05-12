# Copyright (c) 2011 Luca De Feo.

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

cdef extern from "gfp2.c":

  ctypedef struct __mpz_struct:
    pass

  ctypedef struct GF:
    pass

  ctypedef struct GF_params:
    pass

  # Memory management
  void init_GF(GF *x, GF_params* parent)
  void clear_GF(GF *x)
  # Initialization of GF(p,2)
  bint setup_GF(GF_params* field, char* characteristic)
  void free_GF(GF_params*)
  # I/O
  void set_GF(GF* x, char* a, char* b)
  void get_GF(char *a, char *b, GF x)
  # Arithmetic modulo X^2 + 1
  void copy_GF(GF* res, GF x)
  void add_GF(GF *res, GF x, GF y)
  void sub_GF(GF *res, GF x, GF y)
  void neg_GF(GF *res, GF x)
  void scalar_si_GF(GF *res, GF x, long int s)
  void sqr_GF(GF *res, GF x)
  void mul_GF(GF *res, GF x, GF y)
  bint inv_GF(GF *res, GF x)
  bint div_GF(GF *res, GF x, GF y)
  # Miscellaneaous
  int cmp_GF(GF x, GF y)
  bint is_one_GF(GF x)
  bint is_zero_GF(GF x)
  void random_GF(GF *res)
  void print_GF(GF x)

  # Elliptic curve addition
  void mont_ladder(GF *res1x, GF *res1z,
                   GF *res2x, GF *res2z,
                   GF x1, GF z1,
                   GF x2, GF z2,
                   GF dx, GF dz,
                   GF A24)
  # Montgomery point doubling
  void mont_double(GF *resx, GF *resz,
                   GF x, GF z,
                   GF A24)
  # Montgomery point tripling
  void mont_triple(GF *resx, GF *resz,
                   GF x, GF z,
                   GF A24)
  # 3-point ladder to compute P + [t]Q
  # Inputs: t, P, Q, Q - P
  void mont_3ladder(GF* Rx, GF* Rz, __mpz_struct* t,
                    GF Px, GF Pz, GF Qx, GF Pz,
                    GF QPx, GF QPz, GF A24)
  #  Computes [m]P + [n]Q, with P and Q points on the Montgomery curve
  #  with parameters A,B.  Uses Edwards' coordinates for
  #  calculations.
  void shamir(GF *Rx, GF *Ry, GF *Rz,
              GF A, GF B,
              GF Px, GF Py, GF Pz,
              GF Qx, GF Qy, GF Qz,
              __mpz_struct* m, __mpz_struct* n)

  ctypedef struct iso:
    GF u, r

  ctypedef struct iso2:
    pass

  ctypedef struct iso3:
    GF p, X1, Y2, Y1, Y0, coB

  ctypedef struct iso4:
    GF Ap2, iAm2

  # Utility routine to compute (A+2)/4 
  void a24(GF* A24, GF A)
  # Compute an isomorphism of the montgomery curve
  # sending (x,z) to (0,0).
  void isom_comp(iso* iso, GF* iA, GF* iB, GF* iA24,
                 GF A, GF B, GF A24, GF x, GF z)
  # Apply an isomorphism of Montgomery curves
  void isom_apply(GF* X, GF* Y, iso iso, GF x, GF y, GF z)
  # Compute a 2-isogeny of the montgomery curve
  # sending (x,z) to (1,...).
  void iso2_comp(iso2* iso, GF* iA, GF* iB, GF* iA24,
                 GF A, GF B, GF x, GF z)
  # Apply a 2-isogeny of Montgomery curves
  void iso2_apply(GF* X, GF* Y, GF* Z,
                  iso2 iso, GF x, GF y, GF z)
  # Compute a 3-isogeny of the montgomery curve
  void iso3_comp(iso3* iso, GF* iA, GF* iB, GF* iA24,
                 GF A, GF B, GF x, GF z)
  # Apply a 3-isogeny of Montgomery curves
  void iso3_apply(GF* X, GF* Y, GF* Z,
                  iso3 iso, GF x, GF y, GF z)
  # Compute a 4-isogeny of the montgomery curve
  # sending (1,...) to infinity.
  void iso4_comp(iso4* iso, GF* iA, GF* iB, GF A, GF B)
  # Apply a 4-isogeny of Montgomery curves
  void iso4_apply(GF* X, GF* Y, GF* Z, iso4 iso,
                  GF x, GF y, GF z)

  # And finally, the composite isogeny!
  void push_through_iso(GF *A, GF *B, GF *A24,
                        GF Rx, GF Rz,
                        int ell, int *strategy, int h,
                        GF *Px, GF *Py, GF *Pz,
                        GF *Qx, GF *Qy, GF *Qz)
