/*
    Copyright (c) 2011 Luca De Feo.

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#define DEBUG
#define GF_TMP_REGS 9
#define MAX(a,b) (((a)>(b))? (a):(b))


/****************** TYPES *****************/

struct GF_params;
typedef struct GF_params GF_params;

// Elements of GF(p^2)
typedef struct {
  GF_params* parent;
  mpz_t a, b;
} GF;

// The field GF(p^2)
// basically its characteristic and some work registers
struct GF_params {
  mpz_t p, tmp1, tmp2, tmp3;
  GF GFtmp[GF_TMP_REGS];
  gmp_randstate_t state;
  int initialized;
};




/******** IMPLEMENTATION OF GF(p^2) ********/

// Memory management
void init_GF(GF* x, GF_params* parent) {
  mpz_init(x->a);
  mpz_init(x->b);
  x->parent = parent;
}

void clear_GF(GF *x) {
  mpz_clear(x->a);
  mpz_clear(x->b);
}

// Initialization of GF(p,2)
int setup_GF(GF_params* field, const char* characteristic) {
  if (!characteristic) {
    // If p is NULL, use 2^387 * 3^242 - 1 
    // as default value
    mpz_init_set_ui(field->p, 3);
    mpz_pow_ui(field->p, field->p, 242);
    mpz_mul_2exp(field->p, field->p, 387);
    mpz_sub_ui(field->p, field->p, 1);
  } else {
    mpz_init_set_str(field->p, characteristic, 0);
  }
  // Check that the Legendre symbol of -1 is -1
  if (mpz_fdiv_ui(field->p, 4) != 3) {
    mpz_clear(field->p);
    return 0;
  }

  gmp_randinit_default(field->state);
  mpz_init(field->tmp1); mpz_init(field->tmp2); mpz_init(field->tmp3);
  int i;
  for (i = 0 ; i < GF_TMP_REGS ; i++)
    init_GF(&field->GFtmp[i], field);
  field->initialized = 1;

  return 1;
}

void free_GF(GF_params* field) {
  if (field->initialized) {
    int i;
    for (i = 0 ; i < GF_TMP_REGS ; i++)
      clear_GF(&field->GFtmp[i]);
    mpz_clear(field->p); mpz_clear(field->tmp1); mpz_clear(field->tmp2); mpz_clear(field->tmp3);
    field->initialized = 0;
  }
}

// IO of elements
// outputs are strings in base 16
// (to save some bytes)
void set_GF(GF* x, const char* a, const char* b) {
  mpz_set_str(x->a, a, 0);
  mpz_set_str(x->b, b, 0);
}

void get_GF(char *a, char *b, const GF x) {
  gmp_sprintf(a, "%#Zx", x.a);
  gmp_sprintf(b, "%#Zx", x.b);
}

// Arithmetic modulo X^2 + 1

/*
  There seems to be a bug in GMP 4.2.1 that makes mpz_mod give
  unpredictable results when the mpz_t holding the result is the same
  as one of the operands.
*/

void copy_GF(GF* res, const GF x){
  mpz_set(res->a, x.a);
  mpz_set(res->b, x.b);
  res->parent = x.parent;
}

void add_GF(GF *res, const GF x, const GF y) {
  mpz_add(x.parent->tmp1, x.a, y.a);
  mpz_mod(res->a, x.parent->tmp1, x.parent->p);
  mpz_add(x.parent->tmp1, x.b, y.b);
  mpz_mod(res->b, x.parent->tmp1, x.parent->p);
  res->parent = x.parent;
}

void add_GF_ui(GF *res, const GF x, unsigned long int u) {
  mpz_add_ui(x.parent->tmp1, x.b, u);
  mpz_mod(res->b, x.parent->tmp1, x.parent->p);
  mpz_set(res->a, x.a);
  res->parent = x.parent;
}

void sub_GF(GF *res, const GF x, const GF y) {
  mpz_sub(x.parent->tmp1, x.a, y.a);
  mpz_mod(res->a, x.parent->tmp1, x.parent->p);
  mpz_sub(x.parent->tmp1, x.b, y.b);
  mpz_mod(res->b, x.parent->tmp1, x.parent->p);
  res->parent = x.parent;
}

void sub_GF_ui(GF *res, const GF x, unsigned long int u) {
  mpz_sub_ui(x.parent->tmp1, x.b, u);
  mpz_mod(res->b, x.parent->tmp1, x.parent->p);
  mpz_set(res->a, x.a);
  res->parent = x.parent;
}

void neg_GF(GF *res, const GF x) {
  if (mpz_sgn(x.a) == 0)
    mpz_set(res->a, x.a);
  else
    mpz_sub(res->a, x.parent->p, x.a);
  if (mpz_sgn(x.b) == 0)
    mpz_set(res->b, x.b);
  else
    mpz_sub(res->b, x.parent->p, x.b);
  res->parent = x.parent;
}

void scalar_GF(GF *res, const GF x, mpz_t s) {
  mpz_mul(x.parent->tmp1, x.a, s);
  mpz_mod(res->a, x.parent->tmp1, x.parent->p);
  mpz_mul(x.parent->tmp1, x.b, s);
  mpz_mod(res->b, x.parent->tmp1, x.parent->p);
  res->parent = x.parent;
}

void scalar_GF_si(GF *res, const GF x, long int s) {
  mpz_mul_si(x.parent->tmp1, x.a, s);
  mpz_mod(res->a, x.parent->tmp1, x.parent->p);
  mpz_mul_si(x.parent->tmp1, x.b, s);
  mpz_mod(res->b, x.parent->tmp1, x.parent->p);
  res->parent = x.parent;
}

void mul_GF(GF *res, const GF x, const GF y) {
  mpz_add(x.parent->tmp1, x.a, x.b);
  mpz_sub(x.parent->tmp2, y.b, y.a);
  mpz_mul(x.parent->tmp3, x.parent->tmp1, x.parent->tmp2);

  mpz_mul(x.parent->tmp1, x.a, y.b);
  mpz_mul(x.parent->tmp2, y.a, x.b);

  mpz_sub(x.parent->tmp3, x.parent->tmp3, x.parent->tmp1);
  mpz_add(x.parent->tmp3, x.parent->tmp3, x.parent->tmp2);
  mpz_add(x.parent->tmp1, x.parent->tmp1, x.parent->tmp2);

  mpz_mod(res->a, x.parent->tmp1, x.parent->p);
  mpz_mod(res->b, x.parent->tmp3, x.parent->p);
  res->parent = x.parent;
}

void sqr_GF(GF *res, const GF x) {
  mpz_mul(x.parent->tmp1, x.a, x.b);
  mpz_add(x.parent->tmp1, x.parent->tmp1, x.parent->tmp1);

  mpz_add(x.parent->tmp2, x.b, x.a);
  mpz_sub(x.parent->tmp3, x.b, x.a);
  mpz_mul(x.parent->tmp2, x.parent->tmp2, x.parent->tmp3);

  mpz_mod(res->a, x.parent->tmp1, x.parent->p);
  mpz_mod(res->b, x.parent->tmp2, x.parent->p);
  res->parent = x.parent;
}

int inv_GF(GF *res, const GF x) {
  mpz_mul(x.parent->tmp1, x.a, x.a);
  mpz_addmul(x.parent->tmp1, x.b, x.b);
  if (!mpz_invert(x.parent->tmp3, x.parent->tmp1, x.parent->p))
    return 0;

  mpz_mul(x.parent->tmp1, x.b, x.parent->tmp3);
  mpz_neg(x.parent->tmp3, x.parent->tmp3);
  mpz_mul(x.parent->tmp2, x.a, x.parent->tmp3);

  mpz_mod(res->a, x.parent->tmp2, x.parent->p);
  mpz_mod(res->b, x.parent->tmp1, x.parent->p);

  res->parent = x.parent;
  return 1;
}

int div_GF(GF *res, const GF x, const GF y) {
  if (!inv_GF(&x.parent->GFtmp[0], y)) return 0;
  mul_GF(res, x, x.parent->GFtmp[0]);
  
  return 1;
}

// Miscellaneaous
int cmp_GF(const GF x, const GF y) {
  int c = mpz_cmp(x.a, y.a);
  if (c == 0) c = mpz_cmp(x.b, y.b);
  return c;
}

int is_one_GF(const GF x) {
  return (mpz_sgn(x.a) == 0) && (mpz_cmp_ui(x.b, 1) == 0);
}

int is_zero_GF(const GF x) {
  return (mpz_sgn(x.a) == 0) && (mpz_sgn(x.b) == 0);
}

void random_GF(GF *res) {
  mpz_urandomm(res->a, res->parent->state, res->parent->p);
  mpz_urandomm(res->b, res->parent->state, res->parent->p);
}

void print_GF(const GF x) {
  gmp_printf("%Zd*x + %Zd", x.a, x.b);
}




/***** ALGORITHMS RELATED TO ECC **********/

// One step of Montgomery ladder
void mont_ladder(GF *res1x, GF *res1z,
		 GF *res2x, GF *res2z,
		 const GF x1, const GF z1,
		 const GF x2, const GF z2,
		 const GF dx, const GF dz,
		 const GF A24) {
  GF* tmp = x1.parent->GFtmp;

  add_GF(&tmp[4], x1, z1);         // a = (self.x + self.z)
  sub_GF(&tmp[5], x1, z1);         // b = (self.x - self.z)
  sub_GF(&tmp[6], x2, z2);
  add_GF(&tmp[7], x2, z2);
  sqr_GF(&tmp[1], tmp[4]);         // aa = a.square()
  sqr_GF(&tmp[3], tmp[5]);         // bb = b.square()
  sub_GF(&tmp[0], tmp[1], tmp[3]); // e = aa - bb
  mul_GF(&tmp[6], tmp[6], tmp[4]); // da = (P.x - P.z)*a
  mul_GF(&tmp[7], tmp[7], tmp[5]); // cb = (P.x + P.z)*b

  add_GF(&tmp[2], tmp[6], tmp[7]);
  sqr_GF(&tmp[2], tmp[2]);
  mul_GF(&tmp[2], tmp[2], dz);     // x2 = diff.z*(da + cb).square()

  sub_GF(&tmp[8], tmp[6], tmp[7]);
  sqr_GF(&tmp[8], tmp[8]);
  mul_GF(&tmp[8], tmp[8], dx);     // z2 = diff.x*(da - cb).square()

  mul_GF(res1z, A24, tmp[0]);
  add_GF(res1z, *res1z, tmp[3]);
  mul_GF(res1z, *res1z, tmp[0]);   // z1 = e*(bb + self.curve.A24*e))
  mul_GF(res1x, tmp[1], tmp[3]);  // x1 = aa*bb

  copy_GF(res2x, tmp[2]);
  copy_GF(res2z, tmp[8]);
}

/* Montgomery point doubling */
void mont_double(GF *resx, GF *resz,
		 const GF x, const GF z,
		 const GF A24) {
  GF* tmp = x.parent->GFtmp;

  add_GF(&tmp[0], x, z);           // a = (x + z)
  sqr_GF(&tmp[1], tmp[0]);         // aa = a^2
  sub_GF(&tmp[2], x, z);           // b = (x - z)
  sqr_GF(&tmp[3], tmp[2]);         // bb = b^2
  sub_GF(&tmp[4], tmp[1], tmp[3]); // c = aa - bb
  mul_GF(resz, A24, tmp[4]);
  add_GF(resz, *resz, tmp[3]);
  mul_GF(resz, *resz, tmp[4]);   // Z = c (bb + A24 c))
  mul_GF(resx, tmp[1], tmp[3]);  // X = aa bb  
}

/* Montgomery point tripling */
void mont_triple(GF *resx, GF *resz,
		 const GF x, const GF z,
		 const GF A24) {
  GF* tmp = x.parent->GFtmp;

  // Very dirty function, assuming that mont_double uses
  // registers 0 and 2 to store resp. x+z and x-z
  mont_double(&tmp[5], &tmp[6], x, z, A24);

  sub_GF(&tmp[7], tmp[5], tmp[6]);
  add_GF(&tmp[8], tmp[5], tmp[6]);
  mul_GF(&tmp[5], tmp[7], tmp[0]); // da = (x2 - z2)*a
  mul_GF(&tmp[6], tmp[8], tmp[2]); // cb = (x2 + z2)*b

  add_GF(&tmp[7], tmp[5], tmp[6]);
  sqr_GF(&tmp[7], tmp[7]);
  mul_GF(&tmp[7], tmp[7], z);     // X = z*(da + cb)^2

  sub_GF(&tmp[8], tmp[5], tmp[6]);
  sqr_GF(&tmp[8], tmp[8]);
  mul_GF(resz, tmp[8], x);     // Z = x*(da - cb)^2

  copy_GF(resx, tmp[7]);
}

/*
  Converts a Montgomery point to an Edwards point.
  Guarantee: avoids the temporary registers 0-4.
*/
void mont_to_ed(GF* Rx, GF* Ry,
		const GF x, const GF y, const GF z) {
  GF* tmp = x.parent->GFtmp;

  /*
    X = x(x+z) / y(x+z)
    Y = y(x-z) / y(x+z)
  */
  add_GF(&tmp[5], x, z);
  sub_GF(&tmp[6], x, z);
  mul_GF(&tmp[7], y, tmp[5]);
  inv_GF(&tmp[8], tmp[7]);
  mul_GF(&tmp[7], tmp[5], tmp[8]);
  mul_GF(Rx, x, tmp[7]);
  mul_GF(&tmp[7], tmp[6], tmp[8]);
  mul_GF(Ry, y, tmp[7]);
}


// Three-point ladder addition step:
//   P1 = 2 P1
//   P2 = dadd(P1, P2, D2)
//   P3 = dadd(P1, P3, D3)
void mont_tradd(GF *x1, GF *z1,
		GF *x2, GF *z2,
		GF *x3, GF *z3,
		const GF dx2, const GF dz2,
		const GF dx3, const GF dz3,
		const GF A24) {
  GF* tmp = x1->parent->GFtmp;

  // The use of temporary registers in this function is tailored
  // so that it is compatible with
  //       {x1, x2} == {&tmp[0], &tmp[2]}
  //       {z1, z2} == {&tmp[1], &tmp[3]}
  // so that mont_3ladder may safely call it.
  // Be careful when you change indices!

  add_GF(&tmp[4], *x1, *z1);         // a = (self.x + self.z)
  sub_GF(&tmp[5], *x1, *z1);         // b = (self.x - self.z)
  sub_GF(&tmp[6], *x2, *z2);
  add_GF(&tmp[7], *x2, *z2);
  sub_GF(&tmp[0], *x3, *z3);
  add_GF(&tmp[1], *x3, *z3);

  /* P3 */
  mul_GF(&tmp[2], tmp[0], tmp[4]); // da = (P.x - P.z)*a
  mul_GF(&tmp[3], tmp[1], tmp[5]); // cb = (P.x + P.z)*b

  add_GF(&tmp[8], tmp[2], tmp[3]);
  sqr_GF(&tmp[0], tmp[8]);
  mul_GF(x3, tmp[0], dz3);     // x2 = diff.z*(da + cb).square()

  sub_GF(&tmp[8], tmp[2], tmp[3]);
  sqr_GF(&tmp[0], tmp[8]);
  mul_GF(z3, tmp[0], dx3);     // z2 = diff.x*(da - cb).square()

  /* P2 */
  mul_GF(&tmp[6], tmp[6], tmp[4]); // da = (P.x - P.z)*a
  mul_GF(&tmp[7], tmp[7], tmp[5]); // cb = (P.x + P.z)*b

  add_GF(&tmp[2], tmp[6], tmp[7]);
  sqr_GF(&tmp[3], tmp[2]);
  mul_GF(x2, tmp[3], dz2);     // x2 = diff.z*(da + cb).square()

  sub_GF(&tmp[3], tmp[6], tmp[7]);
  sqr_GF(&tmp[8], tmp[3]);
  mul_GF(z2, tmp[8], dx2);     // z2 = diff.x*(da - cb).square()

  /* P1 */
  sqr_GF(&tmp[6], tmp[4]);         // aa = a.square()
  sqr_GF(&tmp[7], tmp[5]);         // bb = b.square()
  sub_GF(&tmp[8], tmp[6], tmp[7]); // e = aa - bb
  mul_GF(&tmp[4], A24, tmp[8]);
  add_GF(&tmp[5], tmp[4], tmp[7]);
  mul_GF(z1, tmp[5], tmp[8]);      // z1 = e*(bb + self.curve.A24*e))
  mul_GF(x1, tmp[6], tmp[7]);      // x1 = aa*bb
}


// 3-point ladder to compute P + [t]Q
// Inputs: t, P, Q, Q - P
void mont_3ladder(GF* Rx, GF* Rz,
		  const mpz_t t,
		  const GF Px, const GF Pz,
		  const GF Qx, const GF Qz,
		  const GF QPx, const GF QPz,
		  const GF A24) {
  GF* tmp = Px.parent->GFtmp;
  // A = 0
  set_GF(&tmp[0], "0", "1"); set_GF(&tmp[1], "0", "0");
  // B = Q
  copy_GF(&tmp[2], Qx); copy_GF(&tmp[3], Qz);
  // C = P
  copy_GF(Rx, Px); copy_GF(Rz, Pz);

  int bit = mpz_sizeinbase(t, 2) - 1;
  for ( ; bit >=0 ; bit--) {
    if (mpz_tstbit(t, bit) == 0) {
      // A = 2A,  B = A + B,  C = A + C
      mont_tradd(&tmp[0], &tmp[1], &tmp[2], &tmp[3], Rx, Rz,
		 Qx, Qz, Px, Pz, A24);
    } else {
      // A = A + B,  B = 2B,  C = B + C
      mont_tradd(&tmp[2], &tmp[3], &tmp[0], &tmp[1], Rx, Rz,
		 Qx, Qz, QPx, QPz, A24);
    }
  }
}

/*
  Computes [m]P + [n]Q, with P and Q points on the Montgomery curve
  with parameters A,B.  Uses Edwards' coordinates for
  calculations.  */
void shamir(GF* Rx, GF* Ry, GF* Rz,
	    const GF A, const GF B,
	    const GF Px, const GF Py, const GF Pz,
	    const GF Qx, const GF Qy, const GF Qz,
	    const mpz_t m, const mpz_t n) {
  
  // some temporary registers
  GF_params* field = A.parent;
  GF* tmp = field->GFtmp;
  // some other dynamically allocated registers 
  GF a, d, aPx, aPy, aQx, aQy, PQx, PQy;

  /* 
     Parameters of the Edwards curve equivalent to this one:
       a = (A+2)/B
       d = (A-2)/B
  */
  inv_GF(&tmp[0], B);
  copy_GF(&tmp[1], A);
  add_GF_ui(&tmp[1], tmp[1], 2);
  init_GF(&a, field);
  mul_GF(&a, tmp[1], tmp[0]);
  sub_GF_ui(&tmp[1], tmp[1], 4);
  init_GF(&d, field);
  mul_GF(&d, tmp[1], tmp[0]);

  /*
    Computing the Ewdards coordinates of P and Q:
      aPx, aPy = Edwards(P)
      aQx, aQy = Edwards(Q)
  */
  init_GF(&aPx, field); init_GF(&aPy, field);
  mont_to_ed(&aPx, &aPy, Px, Py, Pz);
  init_GF(&aQx, field); init_GF(&aQy, field);
  mont_to_ed(&aQx, &aQy, Qx, Qy, Qz);
  
  /*
    Computing P+Q using affine Edwards.
  */
  mul_GF(&tmp[4], aPx, aQx); // tmp4 = C = aPx * aQx
  mul_GF(&tmp[5], aPy, aQy); // tmp5 = D = aPy * aQy
  add_GF(&tmp[0], aPx, aPy); // tmp0 = A = aPx + aPy
  add_GF(&tmp[2], aQx, aQy); // tmp2 = B = aQx + aQy
  mul_GF(&tmp[7], tmp[4], tmp[5]);
  mul_GF(&tmp[6], d, tmp[7]); // tmp6 = E = d * aPx * aQx * aPy * aQy
  sqr_GF(&tmp[8], tmp[6]);
  neg_GF(&tmp[7], tmp[8]);
  add_GF_ui(&tmp[7], tmp[7], 1);
  inv_GF(&tmp[8], tmp[7]); // tmp8 = 1 / (1-E^2)
  add_GF_ui(&tmp[6], tmp[6], 1);
  mul_GF(&tmp[7], a, tmp[4]);
  sub_GF(&tmp[1], tmp[5], tmp[7]);
  mul_GF(&tmp[7], tmp[6], tmp[1]);
  init_GF(&PQy, field);
  mul_GF(&PQy, tmp[7], tmp[8]); // PQy = (1+E)(D - a C) / (1-E^2)
  neg_GF(&tmp[6], tmp[6]);
  add_GF_ui(&tmp[6], tmp[6], 2);
  mul_GF(&tmp[1], tmp[0], tmp[2]);
  sub_GF(&tmp[3], tmp[1], tmp[4]);
  sub_GF(&tmp[1], tmp[3], tmp[5]);
  mul_GF(&tmp[7], tmp[6], tmp[1]);
  init_GF(&PQx, field);
  mul_GF(&PQx, tmp[7], tmp[8]); // PQx = (1-E)(A B - C - D) / (1-E^2)

  int bit = MAX(mpz_sizeinbase(m, 2), mpz_sizeinbase(n, 2)) - 1;
  mpz_set_ui(Rx->a, 0); mpz_set_ui(Ry->a, 0); mpz_set_ui(Rz->a, 0);
  mpz_set_ui(Rx->b, 0); mpz_set_ui(Ry->b, 1); mpz_set_ui(Rz->b, 1);
  Rx->parent = Ry->parent = Rz->parent = Px.parent;
  for ( ; bit >=0 ; bit--){
    /* Double, using projective Edwards */
    add_GF(&tmp[1], *Rx, *Ry);
    sqr_GF(&tmp[0], tmp[1]); // tmp0 = B = (Rx + Ry)^2
    sqr_GF(&tmp[1], *Rx); // tmp1 = C = Rx^2
    sqr_GF(&tmp[2], *Ry); // tmp2 = D = Ry^2
    mul_GF(&tmp[3], a, tmp[1]); // tmp3 = E = a C
    add_GF(&tmp[4], tmp[3], tmp[2]); // tmp4 = F = E + D
    sqr_GF(&tmp[5], *Rz); // tmp5 = H = Rz^2
    scalar_GF_si(&tmp[7], tmp[5], 2);
    sub_GF(&tmp[6], tmp[4], tmp[7]); // tmp6 = J = F - 2H
    sub_GF(&tmp[7], tmp[0], tmp[1]);
    sub_GF(&tmp[8], tmp[7], tmp[2]);
    mul_GF(Rx, tmp[8], tmp[6]); // Rx = (B-C-D) J
    sub_GF(&tmp[7], tmp[3], tmp[2]);
    mul_GF(Ry, tmp[7], tmp[4]); // Ry = (E-D) F
    mul_GF(Rz, tmp[4], tmp[6]); // Rz = F J

    /* Double and Add, using projective Edwards */
    int r = mpz_tstbit(m, bit) | (mpz_tstbit(n, bit) << 1);
    if (r) {
      if (r == 1) {
	mul_GF(&tmp[0], *Rx, aPx); // tmp0 = C = Rx aPx
	mul_GF(&tmp[1], *Ry, aPy); // tmp1 = D = Ry aPy
	add_GF(&tmp[2], aPx, aPy);  // tmp2 = H = aPx + aPy
      } else if (r == 2) {
	mul_GF(&tmp[0], *Rx, aQx); // tmp0 = C = Rx aQx
	mul_GF(&tmp[1], *Ry, aQy); // tmp1 = D = Ry aQy
	add_GF(&tmp[2], aQx, aQy);  // tmp2 = H = aQx + aQy
      } else {
	mul_GF(&tmp[0], *Rx, PQx); // tmp0 = C = Rx PQx
	mul_GF(&tmp[1], *Ry, PQy); // tmp1 = D = Ry PQy
	add_GF(&tmp[2], PQx, PQy);  // tmp2 = H = PQx + PQy
      }
      sqr_GF(&tmp[3], *Rz); // tmp3 = B = Rz^2
      mul_GF(&tmp[5], tmp[0], tmp[1]);
      mul_GF(&tmp[4], d, tmp[5]); // tmp4 = E = d C D
      sub_GF(&tmp[5], tmp[3], tmp[4]); // tmp5 = F = B - E
      add_GF(&tmp[6], tmp[3], tmp[4]); // tmp6 = G = B + E
      add_GF(&tmp[7], *Rx, *Ry);
      mul_GF(&tmp[8], tmp[7], tmp[2]);
      sub_GF(&tmp[7], tmp[8], tmp[0]);
      sub_GF(&tmp[8], tmp[7], tmp[1]);
      mul_GF(&tmp[7], tmp[8], tmp[5]);
      mul_GF(Rx, *Rz, tmp[7]); // Rx = Rz F ((Rx+Ry)H - C - D)
      mul_GF(&tmp[7], a, tmp[0]);
      sub_GF(&tmp[8], tmp[1], tmp[7]);
      mul_GF(&tmp[7], tmp[6], tmp[8]);
      mul_GF(Ry, *Rz, tmp[7]); // Ry = Rz G (D - a C)
      mul_GF(Rz, tmp[5], tmp[6]); // Rz = F G
    }
  }

  /* Convert to Montgomery */
  add_GF(&tmp[0], *Rz, *Ry);
  sub_GF(&tmp[1], *Rz, *Ry);
  mul_GF(Ry, tmp[0], *Rz); // Ry = (Rz+Ry)Rz
  mul_GF(Rz, tmp[1], *Rx); // Rz = (Rz-Ry)Rx
  mul_GF(Rx, tmp[0], *Rx); // Rx = (Rz+Ry)Rx

  clear_GF(&a); clear_GF(&d);
  clear_GF(&aPx); clear_GF(&aPy);
  clear_GF(&aQx); clear_GF(&aQy);
  clear_GF(&PQx); clear_GF(&PQy);
}




/************* ISOGENIES ******************/
typedef struct {
  GF u, r;
} iso;

typedef GF iso2;

typedef struct {
  GF p, p2;
} iso3;

typedef struct {
  GF Ap2;
} iso4;

/* Utility routine to compute (A+2)/4 */
void a24(GF* A24, const GF A) {
  GF_params* field = A.parent;
  GF* tmp = field->GFtmp;

  add_GF_ui(&tmp[0], A, 2);
  mpz_set_ui(field->tmp1, 4);
  mpz_invert(field->tmp2, field->tmp1, field->p);
  scalar_GF(A24, tmp[0], field->tmp2);
}

/*
  Compute an isomorphism of the montgomery curve
  sending (x,z) to (0,0).
*/
void isom_comp(iso* iso, GF* iA, GF* iB, GF* iA24,
	       const GF A, const GF B, const GF A24,
	       const GF x, const GF z) {
  GF* tmp = A.parent->GFtmp;
  
  mont_double(&tmp[1], &tmp[2], x, z, A24);
  div_GF(&tmp[0], tmp[1], tmp[2]);  // P2x = x([2]P) / z([2]P)
  neg_GF(&iso->r, tmp[0]);  // r = -P2x
  scalar_GF_si(&tmp[1], tmp[0], 3);
  add_GF(&tmp[1], tmp[1], A);  // a2 = 3 P2x + A
  mul_GF(&tmp[2], iso->r, z);
  add_GF(&tmp[3], tmp[2], x);
  div_GF(&iso->u, z, tmp[3]);  // u = z / (x - z P2x)
  mul_GF(iA, tmp[1], iso->u);  // iA = a2 u
  mul_GF(iB, B, iso->u);  // iB = B u
  a24(iA24, *iA);
}

/* Apply an isomorphism of Montgomery curves */
void isom_apply(GF* X, GF* Y, GF* unused,
		const iso iso,
		const GF x, const GF y, const GF z) {
  GF* tmp = x.parent->GFtmp;

  mul_GF(&tmp[0], iso.r, z);
  add_GF(&tmp[1], x, tmp[0]);
  if (Y)
    mul_GF(Y, y, iso.u); // Y = y u
  mul_GF(X, tmp[1], iso.u); // X = (x + r z) u
}

/*
  Compute a 2-isogeny of the montgomery curve
  sending (x,z) to (1,...).
*/
void iso2_comp(iso2* iso, GF* iA, GF* iB, GF* iA24,
	       const GF A, const GF B,
	       const GF x, const GF z) {
  GF* tmp = x.parent->GFtmp;

  sub_GF(&tmp[0], x, z);
  sqr_GF(&tmp[1], tmp[0]);
  inv_GF(&tmp[0], tmp[1]);
  mul_GF(&tmp[1], tmp[0], z);
  mul_GF(iso, tmp[1], x); // iA2 = x z / (x-z)^2
  add_GF_ui(&tmp[0], A, 6);
  mul_GF(iB, B, *iso); // iB = B iA2
  mul_GF(iA, tmp[0], *iso); // iA = (A+6) iA2
  a24(iA24, *iA);
}

/* Apply a 2-isogeny of Montgomery curves */
void iso2_apply(GF* X, GF* Y, GF* Z,
		const iso2 iso,
		const GF x, const GF y, const GF z) {
  GF* tmp = x.parent->GFtmp;
  
  sub_GF(&tmp[3], x, z);
  sqr_GF(&tmp[4], tmp[3]);
  if (Y) {
    mul_GF(&tmp[4], x, tmp[4]); // ... X = x iA2 (x - z)^2
    sqr_GF(&tmp[0], x); // Px2 = x^2
    sqr_GF(&tmp[1], z);
    sub_GF(&tmp[2], tmp[0], tmp[1]);
    mul_GF(&tmp[1], y, tmp[2]);
    mul_GF(Y, iso, tmp[1]); // Y = iA2 y (x^2 - z^2)
    mul_GF(Z, z, tmp[0]); // Z = z x^2
  } else {
    mul_GF(Z, z, x); // Z = z x
  }
  mul_GF(X, iso, tmp[4]); // X = iA2 (x - z)^2
}

/*
  Compute a 3-isogeny of the montgomery curve
*/
void iso3_comp(iso3* iso, GF* iA, GF* iB, GF* iA24,
	       const GF A, const GF B,
	       const GF x, const GF z) {
  GF* tmp = x.parent->GFtmp;

  div_GF(&iso->p, x, z);             // p
  sqr_GF(&iso->p2, iso->p);          // p^2
  
  scalar_GF_si(&tmp[3], iso->p, -6); 
  add_GF(&tmp[4], tmp[3], A);
  mul_GF(&tmp[3], tmp[4], iso->p);
  add_GF_ui(&tmp[4], tmp[3], 6);     // (-6p + A)p + 6

  mul_GF(iB, B, iso->p2);      // iB = B p^2
  mul_GF(iA, tmp[4], iso->p);  // iA = ((-6p + A)p + 6)p

  a24(iA24, *iA);
}

/* Apply a 3-isogeny of Montgomery curves */
void iso3_apply(GF* X, GF* Y, GF* Z,
		const iso3 iso,
		const GF x, const GF y, const GF z) {
  GF* tmp = x.parent->GFtmp;

  mul_GF(&tmp[0], z, iso.p);
  sub_GF(&tmp[1], x, tmp[0]); // h = x - p z
                              // if zero, P is in the kernel
  mul_GF(&tmp[2], x, iso.p);
  sub_GF(&tmp[0], tmp[2], z); // rh = x p - z
  sqr_GF(&tmp[3], tmp[0]);
  mul_GF(&tmp[2], x, tmp[3]); // X0 = x (x p - z)^2

  if (Y) {
    mul_GF(&tmp[2], tmp[2], tmp[1]); // X0 *= h
    mul_GF(&tmp[3], x, z);
    sub_GF_ui(&tmp[4], iso.p2, 1);
    mul_GF(&tmp[5], tmp[3], tmp[4]);
    scalar_GF_si(&tmp[3], tmp[5], -2);
    mul_GF(&tmp[4], tmp[0], tmp[1]);
    add_GF(&tmp[5], tmp[3], tmp[4]);
    mul_GF(&tmp[3], tmp[5], tmp[0]); // (rh (rh h + 2xz(1-p^2)))
    sqr_GF(&tmp[7], tmp[1]);
    mul_GF(&tmp[8], tmp[7], z);
    mul_GF(Y, y, tmp[3]); // Y = y (rh (rh h + 2xz(1-p^2)))
    mul_GF(Z, tmp[8], tmp[1]); // Z = h^2 h z
  } else {
    sqr_GF(&tmp[3], tmp[1]);
    mul_GF(Z, tmp[3], z); // Z = h^2 z
  }
  copy_GF(X, tmp[2]);
}

/*
  Compute a 4-isogeny of the Montgomery curve
  sending (1,...) to infinity.
*/
void iso4_comp(iso4* iso, GF* iA, GF* iB, GF* iA24,
	       const GF A, const GF B) {
  GF* tmp = A.parent->GFtmp;

  add_GF_ui(&iso->Ap2, A, 2); // Ap2 = A + 2
  sub_GF_ui(&tmp[0], A, 2);
  neg_GF(&tmp[0], tmp[0]);
  inv_GF(&tmp[2], tmp[0]); // iAm2 = 1 / (2-A)
  add_GF_ui(&tmp[0], A, 6);
  mul_GF(&tmp[1], tmp[0], tmp[2]);
  mul_GF(iB, B, tmp[2]); // iB = B iAm2
  scalar_GF_si(iA, tmp[1], -2); // iA = -2 (A+6) iAm2
  a24(iA24, *iA);
}

/* Apply a 4-isogeny of Montgomery curves */
void iso4_apply(GF* X, GF* Y, GF* Z,
		const iso4 iso,
		const GF x, const GF y, const GF z) {
  GF* tmp = x.parent->GFtmp;

  mul_GF(&tmp[0], x, z); // z1 = x z
  sub_GF(&tmp[2], x, z);
  sqr_GF(&tmp[1], tmp[2]); // x1 = (x - z)^2
  mul_GF(&tmp[2], tmp[0], iso.Ap2); // zA2 = z1 Ap2
  scalar_GF_si(&tmp[3], tmp[0], 4); // fourz = 4 z1
  add_GF(&tmp[6], tmp[1], tmp[2]);
  add_GF(&tmp[4], tmp[1], tmp[3]);
  mul_GF(&tmp[5], tmp[6], tmp[4]); // x0 = (x1+zA2)(x1+fourz)

  if (Y) {
    mul_GF(&tmp[4], x, tmp[1]); // B = x x1
    mul_GF(&tmp[5], tmp[5], tmp[4]); // x0 *= B
    sqr_GF(&tmp[6], z);
    sub_GF(&tmp[7], tmp[0], tmp[6]);
    scalar_GF_si(&tmp[6], tmp[7], 2);
    add_GF(&tmp[7], tmp[6], tmp[1]); // C = x1 + 2(z1 - z^2)
    sqr_GF(&tmp[6], tmp[1]);
    mul_GF(&tmp[8], tmp[2], tmp[3]);
    sub_GF(&tmp[6], tmp[6], tmp[8]); // D = x1^2 - zA2 fourz
    mul_GF(&tmp[8], tmp[7], tmp[6]);
    mul_GF(Y, y, tmp[8]); // Y = y C D
    sqr_GF(&tmp[6], tmp[4]);
    sub_GF_ui(&tmp[7], iso.Ap2, 4);
    neg_GF(&tmp[7], tmp[7]);
    mul_GF(&tmp[8], tmp[6], tmp[7]);
    mul_GF(Z, z, tmp[8]); // Z = z B^2 (4 - Ap2)
  } else {
    sub_GF(&tmp[4], tmp[3], tmp[2]);
    mul_GF(Z, tmp[1], tmp[4]); // Z = x1 (fourz - zA2)
  }
  copy_GF(X, tmp[5]);
}




/******* COMPOSITE ISOGENIES **************/
/* Implementation of a queue */
typedef struct queue_point {
  GF x, z;
  int h;
  struct queue_point *next, *prev;
} queue_point;

#define Q_INIT(q,field) do {	     \
    q = malloc(sizeof(queue_point)); \
    if (q) {			     \
      q->next = q->prev = NULL;	     \
      init_GF(&q->x, field);	     \
      init_GF(&q->z, field);	     \
      q->h = 0;			     \
    }				     \
  } while(0)
#define Q_CLEAR(q) do {   \
      clear_GF(&q->x);	  \
      clear_GF(&q->z);	  \
      free(q);		  \
    } while(0)
#define Q_PUSH(tail,q) do {			\
    tail->next = q;				\
    q->prev = tail;				\
    tail = q;					\
  } while(0)
#define Q_POP(tail,q) do {				\
    q = tail;						\
    tail = tail->prev;					\
    if (tail) {						\
      tail->next = NULL;				\
    }							\
  } while(0)
#define Q_NEXT(q) (q->next)
#define Q_PREV(q) (q->prev)
#define Q_ISHEAD(q) (q->prev==NULL)
#define Q_ISTAIL(q) (q->next==NULL)

// These bits of code are almost identical for 1, 2, 3, 4
// isogenies, thus we "template" them.
#define APPLY_ISOG(apply,obj,lower) do {		\
    for ( tmp = tail ; tmp ; tmp = Q_PREV(tmp)) {	\
      apply(&tmp->x, NULL, &tmp->z, obj,		\
	    tmp->x, tmp->x, tmp->z);			\
      tmp->h = tmp->h - lower;				\
    }							\
    if (Px && Py && Pz)					\
      apply(Px, Py, Pz, obj, *Px, *Py, *Pz);		\
    if (Qx && Qy && Qz)					\
      apply(Qx, Qy, Qz, obj, *Qx, *Qy, *Qz);		\
  } while (0)
#define COMP_ISOG(comp,obj) do {			\
    Q_POP(tail, tmp);					\
    comp(&obj, A, B, A24, *A, *B, tmp->x, tmp->z);	\
    Q_CLEAR(tmp);					\
  } while (0)


union isogenies {
  struct {
    iso d1;
    iso2 d2;
    iso4 d4;
  };
  iso3 d3;
};

/* Push (Px, Py, Pz) and (Qx, Qy, Qz) through the isogeny of kernel
   generated by (Rx, Rz) using the given strategy. */
void push_through_iso(GF *A, GF *B, GF *A24,
		      const GF Rx, const GF Rz,
		      const int ell, int *strategy, int h,
		      GF *Px, GF *Py, GF *Pz,
		      GF *Qx, GF *Qy, GF *Qz) {
  GF_params* field = A->parent;
  int split, i, first = 1;
  union isogenies phi;
  queue_point *tail, *tmp;

  if (ell == 2) {
    init_GF(&phi.d1.u, field);
    init_GF(&phi.d1.r, field);
    init_GF(&phi.d2, field);
    init_GF(&phi.d4.Ap2, field);
  } else {
    init_GF(&phi.d3.p, field);
    init_GF(&phi.d3.p2, field);
  }

  Q_INIT(tail, field);
  copy_GF(&tail->x, Rx);
  copy_GF(&tail->z, Rz);
  tail->h = h;
  while (tail) {
    h = tail->h;
    split = strategy[h];
    // Descend to the floor
    while (h > 1) {
      Q_INIT(tmp, field);
      copy_GF(&tmp->x, tail->x);
      copy_GF(&tmp->z, tail->z);
      for ( i=0 ; i < h - split ; i++) {
	if (ell == 2)
	  mont_double(&tmp->x, &tmp->z,
		      tmp->x, tmp->z, *A24);
	else 
	  mont_triple(&tmp->x, &tmp->z,
		      tmp->x, tmp->z, *A24);
      }
      tmp->h = split;
      Q_PUSH(tail, tmp);
      h = split;
      split = strategy[h];
    }
    // For ell=2, at the first iteration, bring the
    // 2-torsion point to (0,0)
    if (ell == 2 && first) {
      first = 0;
      Q_INIT(tmp, field); // slight abuse
      mont_double(&tmp->x, &tmp->z, tail->x, tail->z, *A24);
      isom_comp(&phi.d1, A, B, A24,
		*A, *B, *A24, tmp->x, tmp->z);
      Q_CLEAR(tmp);
      APPLY_ISOG(isom_apply, phi.d1, 0);
    }
    // Compute and apply the isogeny
    if (ell == 2) {
      COMP_ISOG(iso2_comp, phi.d2);
      APPLY_ISOG(iso2_apply, phi.d2, 1);
    } else {
      COMP_ISOG(iso3_comp, phi.d3);
      APPLY_ISOG(iso3_apply, phi.d3, 1);
    }
  }
  // For ell=2 there is still a 4-isogeny to apply
  if (ell == 2) {
    iso4_comp(&phi.d4, A, B, A24, *A, *B);
    // This works because the queue is empty
    APPLY_ISOG(iso4_apply, phi.d4, 2);
  }

  if (ell == 2) {
    clear_GF(&phi.d1.u);
    clear_GF(&phi.d1.r);
    clear_GF(&phi.d2);
    clear_GF(&phi.d4.Ap2);
  } else {
    clear_GF(&phi.d3.p);
    clear_GF(&phi.d3.p2);
  }
}



/********** FOR TEST PURPOSES *************/
int main(int argc, void ** args) {
  char* characteristic = NULL;
  long tours = 10000;

  if (argc > 1)
    characteristic = args[1];
  if (argc > 2)
    tours = atoi(args[2]);

  GF_params field;
  if (!setup_GF(&field, characteristic)) {
    printf("-1 must be a quadratic non-residue modulo p.\n");
    return -1;
  }
  
  printf("GMP ver. %d.%d.%d\n", __GNU_MP_VERSION, __GNU_MP_VERSION_MINOR, __GNU_MP_VERSION_PATCHLEVEL);

  GF a, b, c, d, e;
  iso2 iso2; iso3 iso3;
  int i;
  init_GF(&a, &field); init_GF(&b, &field); init_GF(&c, &field);
  init_GF(&d, &field); init_GF(&e, &field);
  init_GF(&iso2, &field); init_GF(&iso3.p, &field);
  init_GF(&iso3.p2, &field);
  random_GF(&a); random_GF(&a);
  random_GF(&b); random_GF(&c);
  random_GF(&d); random_GF(&e);
  random_GF(&iso3.p); random_GF(&iso3.p2);

  clock_t time1, time2, ratio;
  time1 = -clock();
  for (i = 0 ; i < tours ; i++) {
    mont_double(&a, &b, c, d, e);
  }
  time1 += clock();
  printf("%d doublings in %f cpu secs.\n", tours, (double)time1/CLOCKS_PER_SEC);

  time2 = -clock();
  for (i = 0 ; i < tours ; i++) {
    iso2_apply(&a, NULL, &b, iso2, c, d, e);
  }
  time2 += clock();
  printf("%d 2-isogenies in %f cpu secs.\n", tours, (double)time2/CLOCKS_PER_SEC);
  
  printf("Ratio: %f\n\n", (double)time1/time2);

  time1 = -clock();
  for (i = 0 ; i < tours ; i++) {
    mont_triple(&a, &b, c, d, e);
  }
  time1 += clock();
  printf("%d triplings in %f cpu secs.\n", tours, (double)time1/CLOCKS_PER_SEC);

  time2 = -clock();
  for (i = 0 ; i < tours ; i++) {
    iso3_apply(&a, NULL, &b, iso3, c, d, e);
  }
  time2 += clock();
  printf("%d 3-isogenies in %f cpu secs.\n", tours, (double)time2/CLOCKS_PER_SEC);

  printf("Ratio: %f\n", (double)time1/time2);
}
