# Copyright (c) 2018 TU Eindhoven
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions
# are met:
#
# 1. Redistributions of source code must retain the above copyright
#    notice, this list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright
#    notice, this list of conditions and the following disclaimer in the
#    documentation and/or other materials provided with the distribution.
#
# 3. Neither the name of the California Institute of Technology nor TU Eindhoven nor
#    the names of its contributors may be used to endorse or promote
#    products derived from this software without specific prior
#    written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
# FOR A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL CALTECH
# OR THE CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
# LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF
# USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
# OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
# SUCH DAMAGE.
r""" Shrinking and expanding sets or in image processing words: erosion and dilation
https://en.wikipedia.org/wiki/Mathematical_morphology
this file gives a rough implementation that gives
 - for dilation: a polytope over approximation of the dilated polytope
 - for erosion: a polytope under approximation of the eroded set
"""
# Created by S. Haesaert 2018/10
from polytope import extreme, qhull
import polytope
import itertools
import numpy as np

def dilate(poly,eps):
    """
    The function dilates a polytope.

    For a given polytope a polytopic over apoproximation of the $eps$-dilated set is computed.
    An e-dilated Pe set of P is defined as:
        Pe = {x+n|x in P ^ n in Ball(e)}
    where Ball(e) is the epsilon neighborhood with norm |n|<e

    The current implementation is quite crude, hyper-boxes are placed over the original vertices
    and the returned polytope is a qhull of these new vertices.
    :param poly: original polytope
    :param eps: positive scalar value with which the polytope is dilated
    :return: polytope
    """

    vertices = extreme(poly)
    dim = len(vertices[0])  # this is the dimensionality of the space
    dil_eps = dim * [[-eps,eps]]
    dil_eps_v = [np.array(n) for n in itertools.product(*dil_eps)]  # vectors with (+- eps,+- eps, +- eps,...)
    new_vertices = []
    for v,d in itertools.product(vertices,dil_eps_v):
        new_vertices += [[v + d]]
        # make box
    return qhull(np.concatenate(new_vertices))


def erode(poly,eps):
    """
    Given a polytope compute eps erosion and give polytope under approximation

    For a given polytope a polytopic under approximation of the $eps$-eroded set is computed.
    An e-eroded Pe set of P is defined as:
        Pe = {x |x+n in P forall n in Ball(e)}
    where Ball(e) is the epsilon neighborhood with norm |n|<e

    The current implementation shifts hyper-planes  with eps over there normal,
    / / / / / / / | +   |
     / / / / / / /| eps |
     / / / / / / /| +   |
     / / / / / / /| eps |
     / / / / / / /| +   |

    :param poly: original polytope
    :param eps: positive scalar value with which the polytope is eroded
    :return: polytope
    """

    A = poly.A
    b = poly.b
    b_e = []
    for A_i,b_i in itertools.product(A,b):
        b_e += [[b_i- eps*np.linalg.norm(A_i,2)]]



    return polytope.Polytope(A,np.array(b_e))
