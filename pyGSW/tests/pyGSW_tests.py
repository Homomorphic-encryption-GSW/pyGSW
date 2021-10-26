from pyGSW.utils import MatrixUtils
from pyGSW.GSW import GSWParams, GSWPublicKey, GSWSecretKey, GSWKeys
from pyGSW.GSW import HomomorphicOperations

import numpy as np

from random import randint
from unittest import main, TestCase

LAMBDA_VALUE = 7  # values 7(3-4 sec per encrypt operation), 8(50-54 sec per encrypt operation) is OK

# GSW paper - https://eprint.iacr.org/2013/340.pdf

class KeyGenTest(TestCase):

    def params_generation(self):

        # Generating GSW parameters
        params = GSWParams.Setup(LAMBDA_VALUE)

        # Generating GSW keys
        keys = GSWKeys(LAMBDA_VALUE)

        # Check that parameter object have all required attributes
        self.assertTrue(hasattr(params.n, "n"))
        self.assertTrue(hasattr(params.q, "q"))
        self.assertTrue(hasattr(params.chi_scale, "chi_scale"))
        self.assertTrue(hasattr(params.m, "m"))
        self.assertTrue(hasattr(params.l, "l"))
        self.assertTrue(hasattr(params.N, "N"))
        self.assertTrue(hasattr(params.L, "L"))

        # Check that key object have all required attributes
        self.assertTrue(hasattr(keys.params.n, "n"))
        self.assertTrue(hasattr(keys.params.q, "q"))
        self.assertTrue(hasattr(keys.params.chi_scale, "chi_scale"))
        self.assertTrue(hasattr(keys.params.m, "m"))
        self.assertTrue(hasattr(keys.params.l, "l"))
        self.assertTrue(hasattr(keys.params.N, "N"))
        self.assertTrue(hasattr(keys.params.L, "L"))

        # Check that with equal security parameter generates equal parameters object
        self.assertTrue(params.n == keys.params.n)
        self.assertTrue(params.q == keys.params.q)
        self.assertTrue(params.chi_scale == keys.params.chi_scale)
        self.assertTrue(params.m == keys.params.m)
        self.assertTrue(params.l == keys.params.l)

    def test_pk_sk_relation(self):
        params = GSWParams.Setup(LAMBDA_VALUE)

        sk = GSWSecretKey.SecretKeyGen(params)
        pk = GSWPublicKey.PublicKeyGen(params, sk)

        # Based on GSW paper, matrix production of secret key vector and public key matrix must be equal to error vector
        check = np.dot(sk.SK, pk.PK) % params.q

        self.assertTrue(np.all(check == (pk.e % params.q)))


class FlattenFuncsTest(TestCase):

    def test_BD_Pof2(self):
        params = GSWParams.Setup(LAMBDA_VALUE)

        vector_a = np.array(np.arange(0, 10))
        vector_b = np.flipud(vector_a)

        # Based on GSW paper, matrix production of BitDecomp(a) result and Powerof2(b) result
        # must be equal to matrix production of a and b vector
        self.assertTrue(np.dot(MatrixUtils.BitDecomp(params, vector_a), MatrixUtils.Powerof2(params, vector_b)) == np.dot(vector_a, vector_b))

    def test_BD_BDI(self):
        params = GSWParams.Setup(LAMBDA_VALUE)

        vector_a = np.array(np.arange(0, 10))

        # Based on GSW paper, BitDecompInverse(BitDecomp(a)) result must be equal to vector a
        self.assertTrue(np.array_equal(MatrixUtils.BitDecompInverse(params, MatrixUtils.BitDecomp(params, vector_a)), vector_a))

    def test_BD_BDI_Pof2_Flatten(self):
        params = GSWParams.Setup(LAMBDA_VALUE)

        vector_a = np.array(np.arange(0, 10))
        vector_b = np.flipud(vector_a)

        vector_a_bit = MatrixUtils.BitDecomp(params, vector_a)

        # Based on GSW paper, for N-size vector a' = dec_to_bit(a) and vector b
        # matrix production of vector a' and Powerof2(b) result must be equal to
        # matrix production of BitDecompInvers(a') result and vector b and must be equal to
        # matrix production of Flatten(a') result and Powerof2(b) result
        first_test = np.dot(vector_a_bit, MatrixUtils.Powerof2(params, vector_b)) == np.dot(
            MatrixUtils.BitDecompInverse(params, vector_a_bit), vector_b)
        second_test = np.dot(MatrixUtils.BitDecompInverse(params, vector_a_bit), vector_b) == np.dot(
            MatrixUtils.Flatten(params, vector_a_bit), MatrixUtils.Powerof2(params, vector_b))

        self.assertTrue(all([first_test, second_test]))


class EncryptionDecryptionTest(TestCase):

    def test_Encryption_and_Decryption(self):
        params = GSWParams.Setup(LAMBDA_VALUE)

        sk = GSWSecretKey.SecretKeyGen(params)
        pk = GSWPublicKey.PublicKeyGen(params, sk)

        # Creating 5 messages in interval [0, 2n-1]
        messages = [randint(0, params.n * 2) for none in range(5)]
        decrypted_messages = []

        # Encrypt every message and trying to decrypt every received ciphertext
        for message in messages:
            ct = pk.Encrypt(params, message)
            decrypted_messages.append(sk.Decrypt(params, ct))

        # Compare messages and decrypted ones
        self.assertTrue(messages == decrypted_messages)


class HomomorphicOperationTest(TestCase):

    def test_Add(self):
        params = GSWParams.Setup(LAMBDA_VALUE)

        test_results = []

        sk = GSWSecretKey.SecretKeyGen(params)
        pk = GSWPublicKey.PublicKeyGen(params, sk)

        # Creating two lists for 5 messages
        messages_a = [randint(1, params.n) for none in range(5)]
        messages_b = [randint(1, params.n) for none in range(5)]

        messages = list(zip(messages_a, messages_b))

        # Add two encrypted messages and trying to decrypt this sum, then compare
        # decrypted sum with plain messages sum
        for message in messages:
            pt_add = message[0] + message[1]

            ct_1 = pk.Encrypt(params, message[0])
            ct_2 = pk.Encrypt(params, message[1])

            ct_add = HomomorphicOperations.Add(params, ct_1, ct_2)

            ct_add_dec = sk.Decrypt(params, ct_add)

            test_results.append(pt_add == ct_add_dec)

        self.assertTrue(all(test_results))

    def test_ConstMult(self):  # low values of consts, mult must be lower 250 for Lambda = 7 !
        params = GSWParams.Setup(LAMBDA_VALUE)

        test_results = []

        # Max value because big noise grow
        MAX_CONST_VALUE = 10

        sk = GSWSecretKey.SecretKeyGen(params)
        pk = GSWPublicKey.PublicKeyGen(params, sk)

        # Creating lists of messages and small constants
        messages = [randint(1, 2 * params.n // MAX_CONST_VALUE) for none in range(5)]
        consts = [randint(1, MAX_CONST_VALUE) for none in range(5)]

        messages = list(zip(messages, consts))

        # Encrypting messages and multiplying received ciphertexts by constants
        # then comparing multiplying results with plain text multiplying
        for message in messages:
            pt_mult = message[0] * message[1]
            ct = pk.Encrypt(params, message[0])

            ct_constmult = HomomorphicOperations.ConstMult(params, ct, message[1])

            ct_mult_dec = sk.Decrypt(params, ct_constmult)

            test_results.append(pt_mult == ct_mult_dec)

        self.assertTrue(all(test_results))


if __name__ == "__main__":
    main()
