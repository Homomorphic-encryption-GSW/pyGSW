""" pyGSW utility functions """

from time import time
from random import randint
from scipy.linalg import block_diag

import numpy as np


start = None
def status(msg):
    """Print time, that past from system work start and input message

    :param msg: printing message
    :type msg: str
    :return: None
    """

    global start
    now = time()
    if start is None:
        start = now
    print("%.4f  %s\n" % (now-start, msg))


class Prime(object):
    """Prime Digits Class"""
    
    @staticmethod
    def is_prime(p):
        """Check whether p is probably prime

        :param p: testing number
        :type p: int
        :return: p prime decision
        :rtype: bool
        """

        for _ in range(16):
            a = randint(1, p-1)
            if pow(a, p-1, p) != 1:
                return False
        return True

    @staticmethod
    def gen_prime(b):
        """Generating a prime p with b bits

        :param b: bit size of generating prime
        :type b: int
        :return: probably prime number with b bits size
        :rtype: int
        """

        p = randint(2**(b-1), 2**b)
        while not Prime.is_prime(p):
            p = randint(2**(b-1), 2**b)
        return p

    @staticmethod
    def generateSafePrime(k):
        """Generating a safe Sophie Germain prime p with k bits

        :param k: bit size of generating Sophie Germain prime
        :type k: int
        :return: probably prime, that is an safe Sophie Germain number with b bits size
        :rtype: int
        """

        p = Prime.gen_prime(k-1)
        sp = 2*p + 1
        while not Prime.is_prime(sp):
            p = Prime.gen_prime(k-1)
            sp = 2*p + 1
        return sp


class MatrixUtils(object):
    """Matrix Flattening functions Class"""
    
    @staticmethod
    def dec_to_bin(x, len):
        """Convert number x in the shape of [LSB ... MSB]

        :param x: converting number
        :type x: int
        :param len: required x bit length
        :type len: int
        :return: l-bit binary convertation of integer x as [LSB ... MSB]
        :rtype: list[bool]
        """

        x = int(x)
        str = bin(x)[2:].zfill(len)
        x = [bool(int(x)) for x in list(str)]
        x = np.asarray(x)
        return np.flipud(x)

    @staticmethod
    def BitDecomp(params, vector):
        """Converting input k-size vector in the shape of (dec_to_bin(v_1)| ... |dec_to_bin(v_k)). Invert

        :param params: GSW scheme parameters
        :type params: GSWParams
        :param vector: decomposable vector
        :type vector: np.array
        :return: k-size vector in the shape of (dec_to_bin(v_1), ..., dec_to_bin(v_k))
        :rtype: np.array
        """

        x = np.zeros(vector.shape[0]*params.l, dtype=int)
        for i in range(0, vector.shape[0]):
            x[(i*params.l):(i*params.l+params.l)] = MatrixUtils.dec_to_bin(vector[i]%params.q, params.l)
        return x

    @staticmethod
    def BitDecompMatrix(params, matrix):
        """Converting input matrix with k-size vectors in the shape of
        (dec_to_bin(v_11)| ... |dec_to_bin(v_1k))
        ...
        (dec_to_bin(v_l1)| ... |dec_to_bin(v_lk))

        :param params: GSW scheme parameters
        :type params: GSWParams
        :param matrix: converting matrix
        :type matrix: np.array
        :return: converted matrix with revers bit representating integers
        :rtype: np.array
        """

        x = []
        for i in range(0, matrix.shape[0]):
            x.append(MatrixUtils.BitDecomp(params, matrix[i]))
        x = np.asarray(x)
        return x

    @staticmethod
    def BitDecompInverse(params, vector):
        """Convert (dec_to_bin(v_1), ..., dec_to_bin(v_k)) input vector to (v_1, ..., v_k)

        :param params: GSW scheme parameters
        :type params: GSWParams
        :param vector: converting vector
        :type vector: np.array
        :return: converted vector with integers
        :rtype: np.array
        """

        x = np.zeros(vector.shape[0]//params.l, dtype=int)
        weight = np.flipud([2**(params.l-z-1) for z in np.arange(params.l)])  # инвертировать???
        for i in range(x.shape[0]):
            x[i] = np.dot(vector[(i*params.l):(i*params.l+params.l)], weight) % params.q    
        return x

    @staticmethod
    def BitDecompInverseMatrix(params, matrix):
        """Convert input matrix
        (dec_to_bin(v_11), ..., dec_to_bin(v_1k))
        ...
        (dec_to_bin(v_l1), ..., dec_to_bin(v_lk))
        to integer matrix
        (v_11, ..., v_1k)
        ...
        (v_l1, ..., v_lk)

        :param params: GSW scheme parameters
        :type params: GSWParams
        :param matrix: bit decomposed matrix
        :type matrix: np.array
        :return: integer matrix
        :rtype: np.array
        """

        x = []
        for i in range(0, matrix.shape[0]):
            x.append(MatrixUtils.BitDecompInverse(params, matrix[i]))
        x = np.asarray(x)
        return x

    @staticmethod
    def Powerof2(params, vector):
        """Convert input k-size vector to k*l-size vector
         (v_1, v_1*2, ..., v_1*2^(l-1), ..., v_k, v_k*2, ..., v_k*2^(l-1))

        :param params: GSW scheme parameters
        :type params: GSWParams
        :param vector: converting integer vector
        :type vector: np.array
        :return: k*l-size vector (v_1, v_1*2, ..., v_1*2^(l-1), ..., v_k, v_k*2, ..., v_k*2^(l-1))
        :rtype: np.array
        """

        weight = np.flipud([2**(params.l-z-1) for z in np.arange(params.l)]) # инвертировать???
        x = []
        for i in range(0, vector.shape[0]):
            x.append(np.dot(vector[i], weight))
        x = np.hstack(np.asarray(x))
        return x

    @staticmethod
    def Flatten(params, vector):
        """Flattening input vector for decreasing vector coefficients

        :param params: GSW scheme parameters
        :type params: GSWParams
        :param vector: flattening vector
        :type vector: np.array
        :return: vector with decreased coefficients
        :rtype: np.array
        """

        return MatrixUtils.BitDecomp(params, MatrixUtils.BitDecompInverse(params, vector))

    @staticmethod
    def FlattenMatrix(params, matrix):
        """Flattening input matrix for decreasing matrix coefficients row by row

        :param params: GSW scheme parameters
        :type params: GSWParams
        :param matrix: flattening matrix
        :type matrix: np.array
        :return: matrix with decreased coefficients
        :rtype: np.array
        """

        return MatrixUtils.BitDecompMatrix(params, MatrixUtils.BitDecompInverseMatrix(params, matrix))

    @staticmethod
    def buildGadget(params):
        """Generating gadget matrix for GSW operations

        :param params: GSW scheme parameters
        :type params: GSWParams
        :return: return matrix (I_n x g), where x - tensor multiplication, g - vector (1, 2, ..., 2^(l-1))
        :rtype: np.array
        """

        g = 2**np.arange(params.l)
        return block_diag(*[g for null in range(params.n)])
