""" pyGSW utility functions """

from time import time
from random import randint
from scipy.linalg import block_diag

import numpy as np


start = None
def status(msg):
    global start
    now = time()
    if start is None:
        start = now
    print("\n%.4f  %s\n" % (now-start, msg))


class Prime(object):
    """ Prime Digits Class """
    
    @staticmethod
    def is_prime(p):
        """ Returns whether p is probably prime """
        for _ in range(16):
            a = randint(1, p-1)
            if pow(a, p-1, p) != 1:
                return False
        return True

    @staticmethod
    def gen_prime(b):
        """ Returns a prime p with b bits """
        p = randint(2**(b-1), 2**b)
        while not Prime.is_prime(p):
            p = randint(2**(b-1), 2**b)
        return p

    @staticmethod
    def generateSafePrime(k):
        """ Return a safe Sophie Germain prime p with k bits """
        p = Prime.gen_prime(k-1)
        sp = 2*p + 1
        while not Prime.is_prime(sp):
            p = Prime.gen_prime(k-1)
            sp = 2*p + 1
        return sp


class MatrixUtils(object):
    """ Matrix Flattening Class """
    
    @staticmethod
    def dec_to_bin(x, len):
        """ Returns l-bit binary representation of integer x as [LSB LSB+1 ... MSB] """
        x = int(x)
        str = bin(x)[2:].zfill(len)
        x = [bool(int(x)) for x in list(str)]
        x = np.asarray(x)
        return np.flipud(x)

    @staticmethod
    def BitDecomp(params, vector):
        x = np.zeros(vector.shape[0]*params.l, dtype=int)
        for i in range(0, vector.shape[0]):
            x[(i*params.l):(i*params.l+params.l)] = MatrixUtils.dec_to_bin(vector[i]%params.q, params.l)
        # print(x)
        return x

    @staticmethod
    def BitDecompMatrix(params, matrix):
        x = []
        for i in range(0, matrix.shape[0]):
            x.append(MatrixUtils.BitDecomp(params, matrix[i]))
        x = np.asarray(x)
        # print(x)
        return x

    @staticmethod
    def BitDecompInverse(params, vector):
        x = np.zeros(vector.shape[0]//params.l, dtype=int)
        weight = np.flipud([2**(params.l-z-1) for z in np.arange(params.l)]) # инвертировать???
        for i in range(x.shape[0]):
            x[i] = np.dot(vector[(i*params.l):(i*params.l+params.l)], weight) % params.q    
        return x

    @staticmethod
    def BitDecompInverseMatrix(params, matrix):
        x = []
        for i in range(0, matrix.shape[0]):
            x.append(MatrixUtils.BitDecompInverse(params, matrix[i]))
        x = np.asarray(x)
        # print(x)
        return x

    @staticmethod
    def Powerof2(params, vector):
        weight = np.flipud([2**(params.l-z-1) for z in np.arange(params.l)]) # инвертировать???
        x = []
        for i in range(0, vector.shape[0]):
            x.append(np.dot(vector[i], weight))
        x = np.hstack(np.asarray(x))
        # print(x)
        return x

    @staticmethod
    def Flatten(params, vector):
        return MatrixUtils.BitDecomp(params, MatrixUtils.BitDecompInverse(params, vector))

    @staticmethod
    def FlattenMatrix(params, matrix):
        return MatrixUtils.BitDecompMatrix(params, MatrixUtils.BitDecompInverseMatrix(params, matrix))

    @staticmethod
    def buildGadget(params):
        # the secret vector [s] is an (n-1)-dimensional vector,
        #   the secret key [t] is -s‖1, an n-dimensional vector
        #
        # the error vector [e] is an m-dimensional vector
        #
        # the matrix [A] is an (n-1)×m matrix (n-1 rows, m = n×l columns)
        #
        # the public key [B] is (   A  ) which is an n×m matrix
        #                       ( sA+e )
        #
        g = 2**np.arange(params.l)
        # print(g)
        return block_diag(*[g for null in range(params.n)])
