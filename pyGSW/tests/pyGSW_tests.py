from pyGSW.utils import MatrixUtils
from pyGSW.pyGSW import GSWParams, GSWPublicKey, GSWSecretKey, GSWKeys
from pyGSW.pyGSW import HomomorphicOperations

import numpy as np

from random import randint
from unittest import main, TestCase

LAMBDA_VALUE = 7  # 7(3-4 sec per encrypt operation), 8(50-54 sec per encrypt operation) is OK


class KeyGenTest(TestCase):

    def params_generation(self):
        params = GSWParams.Setup(LAMBDA_VALUE)

        keys = GSWKeys(LAMBDA_VALUE)

        self.assertTrue(hasattr(params.n, "n"))
        self.assertTrue(hasattr(params.q, "q"))
        self.assertTrue(hasattr(params.chi_scale, "chi_scale"))
        self.assertTrue(hasattr(params.m, "m"))
        self.assertTrue(hasattr(params.l, "l"))
        self.assertTrue(hasattr(params.N, "N"))
        self.assertTrue(hasattr(params.L, "L"))

        self.assertTrue(hasattr(keys.params.n, "n"))
        self.assertTrue(hasattr(keys.params.q, "q"))
        self.assertTrue(hasattr(keys.params.chi_scale, "chi_scale"))
        self.assertTrue(hasattr(keys.params.m, "m"))
        self.assertTrue(hasattr(keys.params.l, "l"))
        self.assertTrue(hasattr(keys.params.N, "N"))
        self.assertTrue(hasattr(keys.params.L, "L"))

        self.assertTrue(params.n == keys.params.n)
        self.assertTrue(params.q == keys.params.q)
        self.assertTrue(params.chi_scale == keys.params.chi_scale)
        self.assertTrue(params.m == keys.params.m)
        self.assertTrue(params.l == keys.params.l)

    def test_pk_sk_relation(self):
        params = GSWParams.Setup(LAMBDA_VALUE)

        sk = GSWSecretKey.SecretKeyGen(params)
        pk = GSWPublicKey.PublicKeyGen(params, sk)

        check = np.dot(sk.SK, pk.PK) % params.q

        self.assertTrue(np.all(check == (pk.e % params.q)))


class FlattenFuncsTest(TestCase):

    def test_BD_Pof2(self):
        params = GSWParams.Setup(LAMBDA_VALUE)

        vector_a = np.array(np.arange(0, 10))
        vector_b = np.flipud(vector_a)

        self.assertTrue(np.dot(MatrixUtils.BitDecomp(params, vector_a), MatrixUtils.Powerof2(params, vector_b)) == np.dot(vector_a, vector_b))

    def test_BD_BDI(self):
        params = GSWParams.Setup(LAMBDA_VALUE)

        vector_a = np.array(np.arange(0, 10))

        self.assertTrue(np.array_equal(MatrixUtils.BitDecompInverse(params, MatrixUtils.BitDecomp(params, vector_a)), vector_a))

    def test_BD_BDI_Pof2_Flatten(self):
        params = GSWParams.Setup(LAMBDA_VALUE)

        vector_a = np.array(np.arange(0, 10))
        vector_b = np.flipud(vector_a)

        vector_a_bit = MatrixUtils.BitDecomp(params, vector_a)

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

        messages = [randint(0, params.n * 2) for none in range(5)]
        decrypted_messages = []

        for message in messages:
            ct = pk.Encrypt(params, message)
            decrypted_messages.append(sk.Decrypt(params, ct))

        self.assertTrue(messages == decrypted_messages)


class HomomorphicOperationTest(TestCase):

    def test_Add(self):
        params = GSWParams.Setup(LAMBDA_VALUE)

        test_results = []

        sk = GSWSecretKey.SecretKeyGen(params)
        pk = GSWPublicKey.PublicKeyGen(params, sk)

        messages_a = [randint(1, params.n) for none in range(5)]
        messages_b = [randint(1, params.n) for none in range(5)]

        messages = list(zip(messages_a, messages_b))

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

        MAX_CONST_VALUE = 10

        sk = GSWSecretKey.SecretKeyGen(params)
        pk = GSWPublicKey.PublicKeyGen(params, sk)

        messages = [randint(1, 2 * params.n // MAX_CONST_VALUE) for none in range(5)]
        consts = [randint(1, MAX_CONST_VALUE) for none in range(5)]

        messages = list(zip(messages, consts))

        for message in messages:
            pt_mult = message[0] * message[1]
            ct = pk.Encrypt(params, message[0])

            ct_constmult = HomomorphicOperations.ConstMult(params, ct, message[1])

            ct_mult_dec = sk.Decrypt(params, ct_constmult)
            print(ct_mult_dec)

            test_results.append(pt_mult == ct_mult_dec)

        self.assertTrue(all(test_results))


if __name__ == "__main__":
    main()
