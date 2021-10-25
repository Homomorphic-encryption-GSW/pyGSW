from pyGSW import GSWKeys
from pyGSW import HomomorphicOperations

from random import randint

LAMBDA_VALUE = 7  # values 7(3-4 sec per encrypt operation), 8(50-54 sec per encrypt operation) is OK


def simple_encryption_decryption():
    """Encryption and decryption example"""

    # To encrypt and decrypt eed to import GSWKeys Class
    # Generating GSW parameters, public key and secret key keypair
    keys = GSWKeys(LAMBDA_VALUE)
    print("GSW parameters: \n", keys.params)
    print("GSW public key: \n", keys.public_key)
    print("GSW secret key: \n\n", keys.secret_key)

    # GSW scheme can encrypt message as integer numbers
    # GSW scheme may work with messages from interval [0, 2^(Lambda+1)] or [0, 2*params.n]
    message = randint(0, pow(2, LAMBDA_VALUE+1))

    # Encrypting generated message
    print("Do some encryption\n")
    print(f"Encrypting message = {message}")
    enc = keys.public_key.Encrypt(keys.params, message)

    # Decrypting operation
    dec = keys.secret_key.Decrypt(keys.params, enc)
    print(f"Decrypted message: {dec}\n")


def simple_homomorphic_operations():
    """Homomorphic operation example"""

    # To perform homomorphic operations need import HomomorphicOperations Class
    # Generating GSW parameters, public key and secret key keypair
    keys = GSWKeys(LAMBDA_VALUE)

    print("Do some add\n")

    # Creating two messages
    # Because sum must be less then 2^(Lambda+1), generate two messages <= 2^(Lambda)
    message_a = randint(1, keys.params.n)
    message_b = randint(1, keys.params.n)

    # Add two encrypted messages and trying to decrypt this sum, then compare
    # decrypted sum with plain messages sum

    # Add two plain messages
    pt_add = message_a + message_b
    print(f"Sum of plain messages: {pt_add}")

    # Encrypting two messages
    ct_a = keys.public_key.Encrypt(keys.params, message_a)
    ct_b = keys.public_key.Encrypt(keys.params, message_b)

    # Performing add operation
    ct_add = HomomorphicOperations.Add(keys.params, ct_a, ct_b)

    # Decrypting sum of ct_a and ct_b
    ct_add_dec = keys.secret_key.Decrypt(keys.params, ct_add)
    print(f"Decrypted sum of encrypted messages: {pt_add}")

    # Other homomorphic operations execute similarly


if __name__ == "__main__":
    simple_encryption_decryption()
    simple_homomorphic_operations()
