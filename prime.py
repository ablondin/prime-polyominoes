r"""
Module to study prime and composed discrete figures

A discrete figure is a connected set of unit squares
without hole. A representation of figure discrete is
given with the help of the Freeman chain code
``F = \{0,1,2,3\}``, where the four letters correspond
to the elementary steps right, up, left, down on the
integral points of the plane.

An homologous morphism ``\varphi`` is a substitution over
``F`` such that ``\varphi(a) = \hat{\varphi(\bar{a})``,
where ``\bar{a}`` is the morphism defined by
``0 \leftrightarrow 2`` and ``1 \leftrightarrow 3``, while
the ``\hat{\cdot}`` operator is defined the composition of
the operator ``\bar{\cdot}`` and the reversal operator
``\tilde{\cdot}``.

One might observe that given an homologous morphism ``\varphi``,
the boundary word ``\varphi(0123)`` is a parallelogram tile.
"""

################################################################################
#                                                                              #
#    Copyright (C) 2013 A. Blondin Masse                                       #
#                                                                              #
#    This program is free software: you can redistribute it and/or modify      #
#    it under the terms of the GNU General Public License as published by      #
#    the Free Software Foundation, either version 3 of the License, or         #
#    (at your option) any later version.                                       #
#                                                                              #
#    This program is distributed in the hope that it will be useful,           #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of            #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             #
#    GNU General Public License for more details.                              #
#                                                                              #
#    You should have received a copy of the GNU General Public License         #
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.     #
#                                                                              #
################################################################################

# Constants
FREEMAN = [0,1,2,3]
HAT     = lambda L: tuple(reversed(map(lambda l: (l + 2) % 4, L)))

# --------- #
# Morphisms #
# --------- #

r"""
Class for representing morphisms
"""
class Morphism(dict):

    def is_homologous(self):
        r"""
        Returns True if self induces a parallelogram tile
        """
        if self is None or len(self) != 4:
            return False
        try:
            BoundaryWord(self[0] + self[1] + self[2] + self[3])
            return True
        except ValueError:
            return False

    def is_one_uniform(self):
        r"""
        Returns True if self 1-uniform, i.e. all images have length
        one
        """
        return all(len(self[letter]) == 1 for letter in self)

    def is_non_trivial_homologous(self):
        r"""
        Returns True if self induces a non trivial parallelogram tile
        """
        return self.is_homologous() and not self.is_one_uniform()

    def copy(self):
        return Morphism(super(Morphism, self).copy())

# -------------- #
# Boundary words #
# -------------- #

r"""
Class for representing boundary words, i.e.
words on the alphabet {0,1,2,3} describin the
boundary of a discrete region.
"""
class BoundaryWord:

    def __init__(self, word):
        r"""
        Creates a boundary word from a list of integers between
        0 and 3 corresponding to the Freeman chain code
        """
        self._word = word
        if not self._is_freeman_word():
            raise ValueError, 'Word must be a list of integers between 0 and 3'
        elif not is_closed(self._word):
            raise ValueError, 'Boundary word is not closed'
        elif not is_simple(self._word[:-1]):
            raise ValueError, 'Boundary word is not simple'

    def _is_freeman_word(self):
        r"""
        Returns True if and only if the boundary word consists
        of letters between 0 and 3
        """
        for letter in self._word:
            if letter not in FREEMAN:
                return False
        return True

    def __repr__(self):
        r"""
        String representation of self
        """
        return 'Boundary word of length %s over {0,1,2,3}' % len(self._word)

    def __getitem__(self, i):
        return self._word[i]

    def turning_number(self):
        r"""
        Returns the turning number of self
        """
        i = n = 0
        while i < len(self._word):
            turn = (self._word[(i + 1) % len(self._word)] - self._word[i]) % 4
            if turn == 1:
                n += 1
            elif turn == 3:
                n -= 1
            i += 1
        return n / 4

    def occurrences_iterator(self, letter, i):
        r"""
        Iterates over all occurrences of the given
        letter starting from index i
        """
        for k in range(i, len(self._word)):
            if self._word[k] == letter:
                yield k

    def is_prime(self):
        r"""
        Returns True if the given boundary word is prime
        """
        return self.factorize() is None

    def factorize(self, verbose=False):
        r"""
        Finds a factorization of the given word
        """
        for i in self.occurrences_iterator(0, 0):
            word = BoundaryWord(self._word[i:] + self._word[:i])
            factorization = word._factorize(Morphism(), 0, verbose)
            if factorization is not None:
                return factorization
        return None

    def _is_non_trivial_factorization(self, morphism):
        r"""
        Returns True if the morphism does not describe the unit square
        nor the whole boundary word
        """
        return morphism.is_non_trivial_homologous() and\
               sum(map(len, morphism.values())) < len(self._word)

    def _factorize(self, morphism, i, verbose=False):
        r"""
        Recursive factorization
        """
        if verbose:
            print morphism, i
        if i >= len(self._word):
            if self._is_non_trivial_factorization(morphism):
                return morphism
            else:
                return None
        else:
            letter = self[i]
            if letter in morphism:
                # We check if we can decode correctly
                k = len(morphism[letter])
                if k > len(self._word) - i or self._word[i:i+k] != morphism[letter]:
                    return None
                else:
                    return self._factorize(morphism.copy(), i + k, verbose)
            else:
                # We try each possible image for morphism[letter]
                for j in self.occurrences_iterator(letter, i):
                    morphism[letter] = tuple(self._word[i:j+1])
                    morphism[(letter + 2) % 4] = HAT(self._word[i:j+1])
                    phi = self._factorize(morphism.copy(), j + 1, verbose)
                    if phi is not None and self._is_non_trivial_factorization(phi):
                        return phi
                return None

# ---------------------- #
# Generating polyominoes #
# ---------------------- #

def is_simple(word):
    r"""
    Returns True if and only if the given word induces a 
    self-avoiding path
    """
    points = {(0,0): True}
    x = y = 0
    for letter in word:
        if letter == 0:
            x += 1
        elif letter == 1:
            y += 1
        elif letter == 2:
            x -= 1
        elif letter == 3:
            y -= 1
        if (x,y) in points:
            return False
        points[(x,y)] = True
    return True 

def is_closed(word):
    r"""
    Returns True if and only if the numbers of 0's and 2'
    are the same and the numbers of 1's and 3's are the
    same
    """
    x = y = 0
    for letter in word:
        if letter == 0:
            x += 1
        elif letter == 1:
            y += 1
        elif letter == 2:
            x -= 1
        elif letter == 3:
            y -= 1
    return x == 0 and y == 0

def simple_words_iterator(length):
    r"""
    Returns an iterator over all simple words of given
    length, i.e.  words inducing self-avoiding paths
    in the discrete plane
    """
    if length == 1:
        for letter in FREEMAN:
            yield (letter,)
    else:
        yielded = {}
        for i in range(1, length):
            first_halves = simple_words_iterator(i)
            for first_half in first_halves:
                second_halves = simple_words_iterator(length - i)
                for second_half in second_halves:
                    word = first_half + second_half
                    if word not in yielded and is_simple(word[:-1]):
                        yielded[word] = True
                        yield word

def boundary_word_iterator(max_length):
    r"""
    Returns an iterator over all boundary words having
    given maximum length
    """
    for length in range(4, max_length + 1, 2):
        for word in simple_words_iterator(length):
            if is_closed(word):
                yield word

# ---- #
# Test #
# ---- #

def is_prime(n):
    from math import sqrt
    for i in range(2, int(sqrt(n)) + 1):
        if n % i == 0:
            return False
    return n >= 2

def test_linear_figures(n):
    for i in range(2, n):
        w = BoundaryWord([0]*i + [1] + [2]*i + [3])
        phi = w.factorize()
        if phi is None:
            if not is_prime(i):
                print i
        else:
            if is_prime(i):
                print i

def prime_polyominoes_iterator(max_length):
    from itertools import ifilter
    return ifilter(lambda w: BoundaryWord(w).is_prime(), boundary_word_iterator(max_length))

def is_lyndon(word):
    found = {}
    for (i,letter) in enumerate(word):
        if letter not in found:
            for j in range(letter):
                if j not in found:
                    return False
            found[letter] = True
            if len(found) == 4:
                return True
    return True

def composed_lyndon_boudary_word_iterator(max_length):
    for w in boundary_word_iterator(max_length):
        bw = BoundaryWord(w)
        if is_lyndon(w) and not bw.is_prime():
            yield w

