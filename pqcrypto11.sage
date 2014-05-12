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

from sage.misc.sage_timeit import sage_timeit
import sage.misc.misc
import paths

load('pqcrypto11.spyx')

################################################################################
#                            EDIT PARAMETERS HERE
################################################################################

# These are some precomputed valid key sizes. The theoreical classical
# security is 1/4 of the bit size (third component of the name). The
# theoretical quantum security is 1/6 of the bit size.
parameters = {
    '2-3-8' : {'lA' : 2, 'lB' : 3, 'eA' : 6, 'eB' : 1, 'f' : 1, 'pm1' : -1},
    '2-3-40' : {'lA' : 2, 'lB' : 3, 'eA' : 22, 'eB' : 15, 'f' : 1, 'pm1' : -1},
    '2-3-256' : {'lA' : 2, 'lB' : 3, 'eA' : 130, 'eB' : 81, 'f' : 22, 'pm1' : -1},
    '2-3-512' : {'lA' : 2, 'lB' : 3, 'eA' : 258, 'eB' : 161, 'f' : 186, 'pm1' : -1},
    '2-3-678' : {'lA' : 2, 'lB' : 3, 'eA' : 341, 'eB' : 218, 'f' : 3, 'pm1' : -1},
    '2-3-768' : {'lA' : 2, 'lB' : 3, 'eA' : 386, 'eB' : 242, 'f' : 2, 'pm1' : -1},
    '2-3-1024' : {'lA' : 2, 'lB' : 3, 'eA' : 514, 'eB' : 323, 'f' : 353, 'pm1' : -1},
    
    '3-5-512' : {'lA' : 3, 'lB' : 5, 'eA' : 161, 'eB' : 110, 'f' : 314, 'pm1' : 1},
    '3-5-512:-1' : {'lA' : 3, 'lB' : 5, 'eA' : 161, 'eB' : 110, 'f' : 736, 'pm1' : -1},
    
    '5-7-32' : {'lA' : 5, 'lB' : 7, 'eA' : 9, 'eB' : 7, 'f' : 16, 'pm1' : -1},
    '5-7-32:-1' : {'lA' : 5, 'lB' : 7, 'eA' : 9, 'eB' : 7, 'f' : 18, 'pm1' : 1},
    '5-7-128' : {'lA' : 5, 'lB' : 7, 'eA' : 55, 'eB' : 46, 'f' : 372, 'pm1' : -1},
    '5-7-512' : {'lA' : 5, 'lB' : 7, 'eA' : 110, 'eB' : 91, 'f' : 284, 'pm1' : -1},
    '5-7-768' : {'lA' : 5, 'lB' : 7, 'eA' : 165, 'eB' : 137, 'f' : 2968, 'pm1' : -1},
    '5-7-1024' : {'lA' : 5, 'lB' : 7, 'eA' : 220, 'eB' : 182, 'f' : 538, 'pm1' : 1},
    
    '11-13-512:+1' : {'lA' : 11, 'lB' : 13, 'eA' : 74, 'eB' : 69, 'f' : 1254, 'pm1' : 1},
    '11-13-512' : {'lA' : 11, 'lB' : 13, 'eA' : 74, 'eB' : 69, 'f' : 384, 'pm1' : -1},
    '11-13-768' : {'lA' : 11, 'lB' : 13, 'eA' : 111, 'eB' : 104, 'f' : 78, 'pm1' : +1},
    '11-13-1024' : {'lA' : 11, 'lB' : 13, 'eA' : 148, 'eB' : 138, 'f' : 942, 'pm1' : +1},
    
    '17-19-512' : {'lA' : 17, 'lB' : 19, 'eA' : 62, 'eB' : 60, 'f' : 120, 'pm1' : -1},
    '17-19-512:+1' : {'lA' : 17, 'lB' : 19, 'eA' : 62, 'eB' : 60, 'f' : 210, 'pm1' : 1},
    '17-19-768' : {'lA' : 17, 'lB' : 19, 'eA' : 94, 'eB' : 90, 'f' : 116, 'pm1' : -1},
    '17-19-1024' : {'lA' : 17, 'lB' : 19, 'eA' : 125, 'eB' : 120, 'f' : 712, 'pm1' : -1},
    
    '23-29-512:-1' : {'lA' : 23, 'lB' : 29, 'eA' : 56, 'eB' : 52, 'f' : 452, 'pm1' : -1},
    '23-29-512' : {'lA' : 23, 'lB' : 29, 'eA' : 56, 'eB' : 52, 'f' : 286, 'pm1' : 1},
    '23-29-768' : {'lA' : 23, 'lB' : 29, 'eA' : 85, 'eB' : 79, 'f' : 132, 'pm1' : -1},
    '23-29-1024' : {'lA' : 23, 'lB' : 29, 'eA' : 113, 'eB' : 105, 'f' : 1004, 'pm1' : -1},
    
    '31-41-512:+1' : {'lA' : 31, 'lB' : 41, 'eA' : 51, 'eB' : 47, 'f' : 1259, 'pm1' : 1},
    '31-41-512' : {'lA' : 31, 'lB' : 41, 'eA' : 51, 'eB' : 47, 'f' : 564, 'pm1' : -1},
    '31-41-768' : {'lA' : 31, 'lB' : 41, 'eA' : 77, 'eB' : 72, 'f' : 166, 'pm1' : +1},
    '31-41-1024' : {'lA' : 31, 'lB' : 41, 'eA' : 103, 'eB' : 95, 'f' : 448, 'pm1' : -1},
    }

# This dictionary is used to tune up the key exchange.  To any prime
# ell is associated a pair giving the relative cost of one
# multiplication by ell and one evaluation of ell-isogeny
# respectively. 
#
# We have some computed some experimental parameters for ell=2,3 using
# an Intel Core Duo U9400. We have left any other prime unoptimized.
#
# To optimize the key exchange for your computer, compile the file
# gfp2.c as follows
#
# > gcc -lgmp gfp2.c
#
# then run it and read the output.
weights = {2 : (2, 1), 3 : (2, 1), 5 : (1, 1), 7 : (1,1),
           11 : (1, 1), 13 : (1, 1), 17 : (1, 1), 19 : (1, 1),
           23 : (1, 1), 29 : (1, 1), 31 : (1, 1), 41 : (1, 1)}


def ss_isogeny(*args, **kw):
    """
    Generate parameters and run a key exchange.

    Input:
      - If no argument is supplied, a default small prime is chosen;
      - If a string is supplied, it is interpreted as the index of one of 
        the predefined primes in the parameters dictionary;
      - Otherwise arguments defining a prime are accepted in the same format
        as for ss_isogeny_gen().

    Optional keywords:
      - verbose (int): if >= 1, print informations and timings about each step
        being executed;
      - measure (bool): passed down to ss_isogeny_exchange().

    Output: the output of ss_isogeny_exchange().
    """
    if 'verbose' in kw:
        old_verbose = misc.get_verbose()
        misc.set_verbose(kw['verbose'])

    if len(args) == 0:
        params = ss_isogeny_gen(**parameters['5-7-32'])
    elif len(args) == 1:
        params = ss_isogeny_gen(**parameters[args[0]])
    else:
        params = ss_isogeny_gen(*args)

    res = ss_isogeny_exchange(*params, measure='measure' in kw and kw['measure'])

    if 'verbose' in kw:
        misc.set_verbose(old_verbose)

    return res
    

################################################################################
#                          SCHEME PARAMETERS (precomputations)
################################################################################
def ss_isogeny_gen(lA, lB, eA, eB, f, pm1):
    """
    Generate public parameters for the given prime.

    Input:
      - lA, lB: prime integers;
      - eA, eB, f: integers;
      - pm1: +1 or -1;
      - lA^eA * lB^eB * f + pm1 must be prime.

    Output:
      tuple (E, lA, eA, PA, QA, lB, eB, PB, QB), where:
      - E is a supersingular curve of cardinality (lA^eA * lB^eB * f)^2;
      - PA and QA are generators of E[lA^eA];
      - PB and QB are generators of E[lB^eB];
      - lA, eA, lB, eB are the same as the inputs.
    """
    prime = "%s^%s * %s^%s" % (lA, eA, lB, eB) + (f != 1)*(" * %s" % f) + (pm1 == 1)*" + 1" + (pm1 == -1)*" - 1"
    p = sage_eval(prime)
    misc.verbose("Using the %s-bits prime p = %s" % (p.nbits(), prime))
    misc.verbose("")

    # Intialize the base field GF(p^2)
    misc.verbose("**** PUBLIC PARAMETERS SET UP ****")
    ct = misc.verbose("Constructing the field GF(p^2)")
    K.<z> = MyGFp2(p)
    P.<X> = K[]
    misc.verbose(t=ct)
    misc.verbose()

    # We compute a supersingular curve using complex multiplication.
    # We need a CM field K such that p is inert in K, if D is the 
    # discriminant of K, this is equivalent to (D/p) = -1.
    # By construction (-1/p) = -1, so the condition above is equivalent
    # to (-D/p) = 1, thus we look for a (small) quadratic residue in GF(p).
    ct = misc.verbose("Looking for a discriminant for which p is inert")
    d = K.base()(-1)
    while not d.is_square():
        d = K.base()(randint(2, 1000))

    D = -d.lift()
    if D % 4 != 1:
        D *= 4
    misc.verbose(t=ct)
    misc.verbose()

    # Now we can compute the Hilbert class polynomial and factor it
    # in GF(p^2)
    ct = misc.verbose("Computing the Hilbert class polynomial for D=%s" % D)
    H = P(hilbert_class_polynomial(D))
    misc.verbose(t=ct)
    misc.verbose()
    ct = misc.verbose("Factoring the Hilbert class polynomial")
    j = H.roots()[-1][0]
    misc.verbose(t=ct)
    misc.verbose()

    # We finally can construct the supersingular elliptic curve
    ct = misc.verbose("Constructing the elliptic curve and taking the right twist")
    E = MontgomeryCurve_from_j(j)
    if ((pm1 == -1 and not (E.A + 2).is_square) or
        (pm1 == 1 and (E.A + 2).is_square)):
        misc.verbose("The curve, has cardinality (p" + (pm1==-1)*"-" + (pm1==1)*"+" + "1)^2, twisting it.")
        E = E.quadratic_twist()
    misc.verbose(t=ct)
    misc.verbose()

    ct = misc.verbose("Checking (probabilistically) that the curve has the expected cardinality")
    for i in range(10):
        a = E.random_point()
        assert((a*(p-pm1)).is_zero())
    misc.verbose(t=ct)
    misc.verbose()

    # We construct the bases of the lA^eA and lB^eB torsion
    def torsion_point(E, cofactor, factor_div_p):
        P = E.zero()
        while (P * factor_div_p).is_zero():
            P = E.random_point().lift() * cofactor
        return P

    def basis(E, cofactor, factor, p):
        P = torsion_point(E, cofactor, factor // p)
        assert((P * factor).is_zero())
        Q = P
        while P.weil_pairing(Q, factor)^(factor // p) == 1:
            Q = torsion_point(E, cofactor, factor // p)
        assert((P * factor).is_zero())
        return (P,Q)

    # Alice's part
    ct = misc.verbose("Constructing Alice's basis")
    PA, QA = basis(E, lB^eB*f, lA^eA, lA)
    misc.verbose(t=ct)
    misc.verbose()

    # Bob's part
    ct = misc.verbose("Constructing Bob's basis")
    PB, QB = basis(E, lA^eA*f, lB^eB, lB)
    misc.verbose(t=ct)
    misc.verbose()

    # We compute the optimal straregies for Alice and Bob
    ct = misc.verbose("Computing Alice's strategy")
    height = eA
    if lA == 2: height -= 2
    strA = paths.optimal_paths(height, *weights[lA], construct=False)
    misc.verbose(t=ct)
    misc.verbose()
    
    ct = misc.verbose("Computing Bob's strategy")
    height = eB
    if lB == 2: height -= 2
    strB = paths.optimal_paths(height, *weights[lB], construct=False)
    misc.verbose(t=ct)
    misc.verbose()

    return E, lA, eA, PA, QA, strA, lB, eB, PB, QB, strB

################################################################################
#                         KEY EXCHANGE
################################################################################
def ss_isogeny_exchange(E, lA, eA, PA, QA, strA, lB, eB, PB, QB, strB, measure=None):
    """
    Perform a key exchange.

    Input:
      - lA, lB: prime integers;
      - eA, eB, f: integers;
      - E a supersingular curve of cardinality (lA^eA * lB^eB * f)^2;
      - PA and QA are generators of E[lA^eA];
      - PB and QB are generators of E[lB^eB];
      - measure: use the sage_timeit module to measure performances
        (default: False)

    Output:
      - The shared key (the j-invariant of a curve isogenous to E);
      - If measure is True, timing information for each phase of the key 
        echange.
    """
    # One run with checks, to see if everything works
    misc.verbose('**** KEY EXCHANGE ****')
    misc.verbose("We run some consistency checks to detect bugs,")
    misc.verbose("so don't take the running times as being accurate.")
    misc.verbose()

    ct = misc.verbose("Randomly generating secret keys.")
    def rand_subgroup(l, e):
        if randrange(0, l+1):
            return (1, randrange(0, l^e))
        else:
            return (l*randrange(0,l^(e-1)), 1)
    mA, nA = rand_subgroup(lA, eA)
    mB, nB = rand_subgroup(lB, eB)

    misc.verbose(t=ct)
    misc.verbose()

    misc.verbose("Generating Alice's public data")
    EA, phiPB, phiQB = keygen_c(PA, QA, mA, nA, lA, strA, PB, QB)
    misc.verbose(t=ct)
    misc.verbose()
    ct = misc.verbose("Generating Bob's public data")
    EB, phiPA, phiQA = keygen_c(PB, QB, mB, nB, lB, strB, PA, QA)
    misc.verbose(t=ct)
    misc.verbose()

    ct = misc.verbose("Computing shared key on Alice's side")
    EsA, _, _ = keygen_c(phiPA, phiQA, mA, nA, lA, strA)
    misc.verbose(t=ct)
    misc.verbose()
    ct = misc.verbose("Computing shared key on Bob's side")
    EsB, _, _ = keygen_c(phiPB, phiQB, mB, nB, lB, strB)
    misc.verbose(t=ct)
    misc.verbose()

    if EsA.j_invariant() != EsB.j_invariant():
        raise RuntimeError, "ERROR: the shared keys don't match! Here's the secret keys:\n\tmA = %d\n\tnA = %d\n\tmB = %d\n\tnB = %d\n" % (mA,nA,mB,nB)

    # Now we measure performances
    timings = None
    if measure is not None and measure:
        misc.verbose('**** TIMINGS ****')
        misc.verbose("Now we measure the real performances.")
        misc.verbose("(This may take some time)")
        misc.verbose()

        context = globals()
        context.update(locals())

        if measure is True:
            repeat = 3
        else:
            try:
                repeat = int(measure)
            except:
                repeat = 3

        misc.verbose("Alice round 1")
        A1 = sage_timeit('keygen_c(PA, QA, mA, nA, lA, strA[:], PB, QB, checks=False)', context, repeat=repeat)
        misc.verbose(A1)
        misc.verbose("Alice round 2")
        A2 = sage_timeit('keygen_c(phiPA, phiQA, mA, nA, lA, strA[:], checks=False)', context, repeat=repeat)
        misc.verbose(A2)
        misc.verbose()

        misc.verbose("Bob round 1")
        B1 = sage_timeit('keygen_c(PB, QB, mB, nB, lB, strB[:], PA, QA, checks=False)', context, repeat=repeat)
        misc.verbose(B1)
        misc.verbose("Bob round 2")
        B2 = sage_timeit('keygen_c(phiPB, phiQB, mB, nB, lB, strB[:], checks=False)', context, repeat=repeat)
        misc.verbose(B2)
        misc.verbose()
        
        timings = A1,A2,B1,B2
        
    return EsA.j_invariant(), timings
