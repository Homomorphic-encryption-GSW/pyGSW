# from pyGSW import pyGSW
from pyGSW.utils import Prime, MatrixUtils, status

import math
import numpy as np


STANDARD_SECURITY_VALUE = 8


class GSWParams(object):
    """Contains GSW scheme parameters

    :arg n: responsible for keys size as 2^Lambda
    :type n: int
    :arg q: responsible for security of scheme as a module of LWE ring
    :type q: int
    :arg chi_scale: error vector distribution. Must be small
    :type chi_scale: float
    :arg m: responsible for public key size
    :type m: int
    :arg l: max bit length of q
    :type l: int
    :arg N: standard GSW parameter, duplicate m
    :type N: int
    :arg L: Deep of homomorphic operations
    :type L: int
    :return: None
    """

    def __init__(self, n: int, q: int, chi_scale, m: int, l: int, N: int, L: int):
        self.n = n
        self.q = q
        self.chi_scale = chi_scale  # q/8 relation???
        self.m = m
        self.l = l
        self.N = N
        self.L = L

    def __str__(self) -> str:
        return "n: {}, q: {}, chi_scale: {}, m: {}, l: {}, N: {}, L: {}".format(
            self.n, self.q, self.chi_scale, self.m, self.l, self.N, self.L)

    @staticmethod
    def Setup(Lambda: int, L: int=10):
        """Generate and setup GSW scheme parameters

        :param Lambda: GSW security parameter. Allows encrypt messages <= 2^(Lambda + 1)
        :type Lambda: (int)
        :param L: Deep of homomorphic operations
        :type L: (int)
        :return: GSW scheme parameters for given security parameter Lambda.
        :rtype: (GSWParams)
        """

        status("Setup GSW parameters")

        n = pow(2, Lambda)

        q = Prime.generateSafePrime(2*Lambda)

        chi_scale = 1.0  # must be around 1.0

        m = n * (math.floor(math.log(q, 2)) + 1)
        l = math.ceil(math.log(q, 2))
        N = n * l

        return GSWParams(n, q, chi_scale, m, l, N, L)


class GSWSecretKey(object):
    """Contains GSW Secret key

    :arg s: secret key in the shape of (t | 1)
    :type s: np.array
    :arg t: secret vector with random value from Z_q
    :type t: np.array
    :arg v: evaluation vector in the shape of (s_1, s_1*2, ..., s_1*2^(l-1), ..., s_k, s_k*2, ..., s_k*2^(l-1))
    :type v: np.array
    :return: None
    """

    def __init__(self, s, t, v):
        self.SK = s
        self.t = t
        self.v = v

    def __str__(self) -> str:
        return f"SK = {self.SK}\n \
                \rt = {self.t}\n \
                \rv = {self.v}"

    @staticmethod
    def SecretKeyGen(params):
        """Generating GSW secret key from input parameters

        :param params: GSW scheme parameters
        :type params: GSWParams
        :return: GSW secret key of GSW parameters
        :rtype: GSWSecretKey
        """

        status("Generating GSW Secret key")

        t = np.random.randint(0, params.q, params.n-1, dtype=np.int64)
        s = np.hstack((t, (np.array([1]))))
        v = MatrixUtils.Powerof2(params, s) % params.q

        return GSWSecretKey(s, t, v)

    def Decrypt(self, params, ciphertext):
        """Decrypting input ciphertext with bonded GSW secret key

        :param params: GSW scheme parameters
        :type params: GSWParams
        :param ciphertext:
        :type ciphertext: np.array
        :return: ciphertext decrypted value as best distance number
        :rtype: int
        """

        status("Decrypting message")

        message = np.dot(self.SK, ciphertext) % params.q

        sg = MatrixUtils.Powerof2(params, self.SK) % params.q

        div = np.rint((message / sg).astype(np.float)).astype(np.int64)

        modes = np.unique(div, return_counts=True)

        modes = sorted(zip(modes[0], modes[1]), key=lambda t: -t[1])
        best_number = 0
        best_dist = float('inf')
        for mu, _ in modes:
            dist = (message - mu*sg) % params.q
            dist = np.minimum(dist, params.q - dist)

            dist = np.dot(dist, dist)
            if dist < best_dist:
                best_number = mu
                best_dist = dist

        return best_number


class GSWPublicKey(object):
    """Contains GSW Public key

    :arg A: GSW public key as a big matrix
    :type A: np.array
    :arg e: GSW error vector. Values in must be small, because after homomorphic operations error little grow
    :type e: np.array
    :return: None

    """

    def __init__(self, A, e):
        self.PK = A
        self.e = e
    
    def __str__(self) -> str:
        return f"PK = {self.PK}\n \
                \re = {self.e}"

    @staticmethod
    def PublicKeyGen(params, sk):
        """Generating GSW public key from input parameters and secret key

        :param params: GSW scheme parameters
        :type params: GSWParams
        :param sk: GSW secret key of GSW parameters
        :type sk: np.array
        :return: GSW public key of GSW parameters. Bonded with input secret key
        :rtype: GSWPublicKey
        """

        status("Generating GSW Public key")

        B = np.random.randint(0, params.q, (params.n-1, params.m), dtype=np.int64)
        e = np.rint(np.random.normal(scale=params.chi_scale, size=params.m)).astype(np.int)

        b = np.add(np.dot(sk.t, B), e) % params.q

        A = np.vstack((-B, b)) % params.q  # B = (n-1) x m; b = 1 x m; A = n x m

        return GSWPublicKey(A, e)

    def Encrypt(self, params, message):
        """Encrypting input integer message with bonded GSW public key

        :param params: GSW scheme parameters
        :type params: GSWParams
        :param message: encrypting integer message. Must be <= 2*(Lambda+1)
        :type message: int
        :return: ciphertext matrix
        :rtype: np.array
        """

        status("Encrypting message")

        R = np.random.randint(2, size=(params.m, params.m), dtype=np.int64)
        G = MatrixUtils.buildGadget(params)

        C = (message*G + np.dot(self.PK, R)) % params.q

        return C


class GSWKeys(object):
    """Construction, containing GSW parameters and bonded secret and public keys of security parameter Lambda

    :arg Lambda: GSW security parameter. Allows encrypt messages <= 2^(Lambda + 1)
    :return: None
    """

    def __init__(self, Lambda):
        self.params = GSWParams.Setup(Lambda)
        self.secret_key = GSWSecretKey.SecretKeyGen(self.params)
        self.public_key = GSWPublicKey.PublicKeyGen(self.params, self.secret_key)


class HomomorphicOperations(object):
    """Contain follow homomorphic operations over ciphertext matrices: Add, Constant Multiplication"""

    @staticmethod
    def Add(params, ciphertext_1, ciphertext_2):
        """Sum two equal shapes ciphertexts matrix, encrypted with equal key

        :param params: GSW scheme parameters
        :type params: GSWParams
        :param ciphertext_1: first ciphertext matrix
        :type ciphertext_1: np.array
        :param ciphertext_2: second ciphertext
        :type ciphertext_2: np.array
        :return: Sum of two input ciphertexts matrix
        :rtype: np.array
        """

        ct1_plus_ct2 = (ciphertext_1 + ciphertext_2) % params.q

        return ct1_plus_ct2

    @staticmethod
    def ConstMult(params, ciphertext, const):
        """Multiplication of input constant and ciphertext

        :param params: GSW scheme parameters
        :type params: GSWParams
        :param ciphertext: multiplied ciphertext matrix
        :type ciphertext: np.array
        :param const: multiplied constant
        :type const: int
        :return: new ciphertext as (const x ciphertext)
        :rtype: np.array
        """

        ct_x_const = np.copy(ciphertext)
        for _ in range(const - 1):
            ct_x_const = (ct_x_const + ciphertext) % params.q

        return ct_x_const

    @staticmethod
    def Mult(params, ciphertext_1, ciphertext_2):  # DO NOT WORK AT ALL! :(
        """Should multiply two input ciphertext, but...

        :param params: GSW scheme parameters
        :type params: GSWParams
        :param ciphertext_1: first ciphertext matrix
        :type ciphertext_1: np.array
        :param ciphertext_2: second ciphertext matrix
        :type ciphertext_2: np.array
        :return: production of two input ciphertext matrix
        :rtype: np.array
        """

        gadget = MatrixUtils.buildGadget(params.l, params.n)
        matrix = MatrixUtils.BitDecompMatrix(params, ciphertext_2.T)

        ca_x_cb = np.dot(ciphertext_1, matrix) % params.q

        return ca_x_cb


if __name__ == "__main__":
    """Encryption and decryption example"""

    keys = GSWKeys()
    print("GSW parameters: \n", keys.params)
    print("GSW public key: \n", keys.public_key)
    print("GSW secret key: \n", keys.secret_key)

    print("Do some encryption")
    print("Encrypting message = 100")
    enc = keys.public_key.Encrypt(keys.params, 100)
    print("Decrypted message:")
    print(keys.secret_key.Decrypt(keys.params, enc))
