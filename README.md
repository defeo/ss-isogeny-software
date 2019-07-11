Quantum-Resistant Cryptosystems from Supersingular Elliptic Curve Isogenies
===========================================================================

Copyright 2011-2016 Luca De Feo <http://defeo.lu/>.

This software implements the cryptosystem described in

> D. Jao and L. De Feo, Towards quantum-resistant cryptosystems from
> supersingular elliptic curve isogenies. *Post-Quantum Cryptography*,
> Nov 2011, Taipei, Taiwan. Springer, LNCS 7071, pp. 19-34, 2011.

> L. De Feo, D. Jao and J. Plût, Towards quantum-resistant
> cryptosystems from supersingular elliptic curve isogenies.  *Journal
> of Mathematical Cryptology*, 8(3), pp. 209-247. De Gruyter, 2014.


**WARNING:** This code is obsolete. For a modern treatment, please see
the official code for the [NIST candidate SIKE](https://sike.org/),
and the additional implementations referenced
[here](https://sike.org/#implementation).


Installation
------------

Just clone or
[download](https://github.com/defeo/ss-isogeny-cryptosystem/archive/master.zip)
this repo.

You will need a recent version of [Sage](http://sagemath.org/) and a C
compiler. This version has been tested with Sage 6.10 and gcc 5.2.1.


Usage
-----

In a Sage shell type

	sage: load('pqcrypto11.sage')

Some predefined key sizes are stored in a string-indexed dictionary
called 'parameters'. Read `pqcrypto11.sage` to find them out.

Public data for a cryptosystem are generated via a call to
`ss_isogeny_gen`. For example, to obtain parameters relative to a
40-bit prime, type

	sage: set_verbose(1)
	sage: pdata = ss_isogeny_gen(**parameters['2-3-40'])

The key exchange is performed by `ss_isogeny_exchange`. Type

	sage: ss_isogeny_exchange(*pdata)
	sage: set_verbose(0)

The function `ss_isogeny` runs both previous functions in one. The
previous sequence of commands is equivalent to

	sage: ss_isogeny('2-3-40', verbose=1)

Additional parameters can be passed to these functions, read
`pqcrypto11.sage`.

**NOTE:** The file `gfp2.c` can be compiled as a standalone program
with

	gcc -lgmp gfp2.c

Then it can be run to gather estimates on the running times of
doublings, triplings, 2 and 3-isogeny evaluations. These data can be
used to tune up (via the dictionary "weights" in `pqcrypto11.sage`)
the key exchange algorithm.


Thanks
------

Many thanks to those who have helped in testing and fixing this
software.

- David Jao,
- Jérôme Plût,
- Erik Nellessen.
- Adarsh Saraf,
- Srinath,
- Miha Marolt @miham
