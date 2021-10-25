# from pyGSW.pyGSW import 
# from utils import MatrixUtils #### ?
# from pyGSW.pyGSW.pyGSW import GSWParams, GSWPublicKey, GSWSecretKey, GSWKeys
# from pyGSW import HomomorphicOperations

from pyGSW.pyGSW.utils import *

from pyGSW.pyGSW.pyGSW import GSWParams, GSWPublicKey, GSWSecretKey, GSWKeys
from pyGSW.pyGSW.pyGSW import HomomorphicOperations

import numpy as np

from random import randint
from unittest import main, TestCase


def test_BD_Pof2(params):
    vector_a = np.array(np.arange(0,10))
    vector_b = np.flipud(vector_a)

    print(f"Test of BitDecomp and Powerof2 funcs: \
    {np.dot(MatrixUtils.BitDecomp(params, vector_a), MatrixUtils.Powerof2(params, vector_b)) == np.dot(vector_a, vector_b)}")


def test_BD_BDI(params):
    vector_a = np.array(np.arange(0,10))

    print(f"Test of BitDecomp and BitDecompInverse funcs: \
    {np.array_equal(MatrixUtils.BitDecompInverse(params, MatrixUtils.BitDecomp(params, vector_a)), vector_a)}")


def test_BD_BDI_Pof2_Flatten(params):
    vector_a = np.array(np.arange(0,10))
    vector_b = np.flipud(vector_a)

    vector_a_bit = MatrixUtils.BitDecomp(params, vector_a)

    first_test = np.dot(vector_a_bit, MatrixUtils.Powerof2(params, vector_b)) == np.dot(MatrixUtils.BitDecompInverse(params, vector_a_bit), vector_b)
    second_test = np.dot(MatrixUtils.BitDecompInverse(params, vector_a_bit), vector_b) == np.dot(MatrixUtils.Flatten(params, vector_a_bit), MatrixUtils.Powerof2(params, vector_b))

    print(f"Test of BitDecomp, BitDecompInverse, Powerof2 and Flatten funcs: \
    {all([first_test, second_test])}")


def test_KeyGen(params):
    SK = GSWSecretKey.SecretKeyGen(params)
    PK = GSWPublicKey.PublicKeyGen(params, SK)

    check = np.dot(SK.SK, PK.PK) % params.q
    okay = np.all(check == (PK.e % params.q))
    print(f"Test of KeyGen funcs: {okay}")


def test_Encryption_and_Decryption(params):
    SK = GSWSecretKey.SecretKeyGen(params)
    PK = GSWPublicKey.PublicKeyGen(params, SK)

    messages = [randint(0, params.n*2) for none in range(5)]
    decrypted_messages = []

    for message in messages:
        print(message)
        ct = GSWPublicKey.Encrypt(params, message)
        decrypted_messages.append(GSWSecretKey.Decrypt(params, ct))
        print(decrypted_messages)

    okay = messages == decrypted_messages

    print(f"Test of Encryption and Decryption funcs: {okay}")


def test_Add(params):
    test_results = []

    SK = GSWSecretKey.SecretKeyGen(params)
    PK = GSWPublicKey.PublicKeyGen(params, SK)

    messages_a = [randint(1, params.n) for none in range(5)]
    messages_b = [randint(1, params.n) for none in range(5)]

    messages = list(zip(messages_a, messages_b))
    print(messages)

    for message in messages:
        pt_add = message[0] + message[1]

        ct_1 = GSWPublicKey.Encrypt(params, message[0])
        ct_2 = GSWPublicKey.Encrypt(params, message[1])

        ct_add = HomomorphicOperations.Add(params, ct_1, ct_2)

        ct_add_dec = GSWSecretKey.Decrypt(params, ct_add)

        test_results.append(pt_add == ct_add_dec)

    print(f"Test of Add func: {all(test_results)}")


def test_ConstMult(params): # low values of consts, mult must be lower 250 for Lambda = 7!
    test_results = []

    MAX_CONST_VALUE = 10

    SK = GSWSecretKey.SecretKeyGen(params)
    PK = GSWPublicKey.PublicKeyGen(params, SK)

    messages = [randint(1, 2*params.n//MAX_CONST_VALUE) for none in range(5)]
    consts = [randint(1, MAX_CONST_VALUE) for none in range(5)]

    messages = list(zip(messages, consts))
    print(messages)

    for message in messages:
        pt_mult = message[0] * message[1]
        print(message[0], message[1])
        ct = GSWPublicKey.Encrypt(params, PK, message[0])

        ct_constmult = HomomorphicOperations.ConstMult(params, ct, message[1])

        ct_mult_dec = GSWSecretKey.Decrypt(params, ct_constmult)
        print(ct_mult_dec)

        test_results.append(pt_mult == ct_mult_dec)

    print(test_results)
    print(f"Test of Add func: {all(test_results)}")


def test():
    LAMBDA_VALUE = 7 ### 7(3-4 sec per ecrypt operation), 8(50-54 sec per encrypt operation) is OK
    
    params = GSWParams.Setup(LAMBDA_VALUE) 
    print(params)
    
    # test_BD_Pof2(params)
    # test_BD_BDI(params)
    # test_BD_BDI_Pof2_Flatten(params)
    # test_KeyGen(params)
    # test_Encryption_and_Decryption(params)
    test_Add(params)
    # test_ConstMult(params)

    vector = np.array([100,2,3,4,5])
    vector = np.array([[1,2,3,4,5],[1,2,3,4,5]])

