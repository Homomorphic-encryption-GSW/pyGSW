# from pyGSW import pyGSW
from utils import Prime, MatrixUtils, status

import math
import numpy as np


STANDARD_SECURITY_VALUE = 8


class GSWParams(object):
    """  """
    
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
        """  """

        status("Setup GSW parameters")

        n = pow(2, Lambda)

        q = Prime.generateSafePrime(2*Lambda)

        chi_scale = 1.0 # must be around 1.0

        m = n * (math.floor(math.log(q, 2)) + 1)
        l = math.ceil(math.log(q, 2))
        N = n * l

        return GSWParams(n, q, chi_scale, m, l, N, L)


class GSWSecretKey(object):
    """  """
    
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
        """  """

        status("Generating GSW Secret key")

        t = np.random.randint(0, params.q, params.n-1, dtype=np.int64)
        s = np.hstack((t, (np.array([1]))))
        v = MatrixUtils.Powerof2(params, s) % params.q

        return GSWSecretKey(s, t, v)

    def Decrypt(self, params, ciphertext):
        """  """

        status("Decrypting message")

        message  = np.dot(self.SK, ciphertext) % params.q 
        gadget = MatrixUtils.buildGadget(params)

        # sg = np.dot(self.SK, gadget) % params.q # Powerof2(SK)
        sg = MatrixUtils.Powerof2(params, self.SK) % params.q
        # print(np.array_equal(sg, sg_2))

        div = np.rint((message / sg).astype(np.float)).astype(np.int64)

        modes = np.unique(div, return_counts=True)

        modes = sorted(zip(modes[0], modes[1]), key = lambda t: -t[1])
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
    """  """
    
    def __init__(self, A, e):
        self.PK = A
        self.e = e
    
    def __str__(self) -> str:
        return f"PK = {self.PK}\n \
                \re = {self.e}"

    @staticmethod
    def PublicKeyGen(params, SK):
        """  """

        status("Generating GSW Public key")

        B = np.random.randint(0, params.q, (params.n-1, params.m), dtype=np.int64)
        e = np.rint(np.random.normal(scale=params.chi_scale, size=params.m)).astype(np.int)

        b = np.add(np.dot(SK.t,B), e) % params.q

        A = np.vstack((-B, b)) % params.q # B = (n-1) x m; A = n x m

        return GSWPublicKey(A, e)

    def Encrypt(self, params, message):
        """  """

        status("Encrypting message")

        R = np.random.randint(2, size=(params.m, params.m), dtype=np.int64)
        G = MatrixUtils.buildGadget(params)

        C = (message*G + np.dot(self.PK, R)) % params.q

        return C


class GSWKeys(object):
    """  """

    def __init__(self):
        self.params = GSWParams.Setup(STANDARD_SECURITY_VALUE)
        self.secret_key = GSWSecretKey.SecretKeyGen(self.params)
        self.public_key = GSWPublicKey.PublicKeyGen(self.params, self.secret_key)


class HomomorphicOperations(object):
    """  """

    @staticmethod
    def Add(params, ciphertext_1, ciphertext_2):
        """  """

        ct1_plus_ct2 = (ciphertext_1 + ciphertext_2) % params.q

        return ct1_plus_ct2

    @staticmethod
    def ConstMult(params, ciphertext, const):
        """  """

        ct_x_const = np.copy(ciphertext)
        for _ in range(const - 1):
            ct_x_const = (ct_x_const + ciphertext) % params.q

        return ct_x_const

    @staticmethod
    def Mult(params, ciphertext_1, ciphertext_2): # DO NOT WORK AT ALL! :(
        """  """

        gadget = MatrixUtils.buildGadget(params.l, params.n)
        matrix = MatrixUtils.BitDecompMatrix(params, ciphertext_2.T)

        new_mat = np.dot(gadget, matrix.T)
        print(f"new_mat == cb: {np.array_equal(new_mat, ciphertext_2)}")

        ca_x_cb = np.dot(ciphertext_1, matrix) % params.q

        return ca_x_cb



if __name__ == "__main__":
    keys = GSWKeys()
    print(keys.params)
    print(keys.public_key)
    print(keys.secret_key)

    enc = keys.public_key.Encrypt(keys.params, 100)
    print(enc)
    print(keys.secret_key.Decrypt(keys.params, enc))

    # python pyGSW\pyGSW\pyGSW.py 
