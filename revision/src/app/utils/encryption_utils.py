import base64

from cryptography.hazmat.backends import default_backend
from cryptography.hazmat.primitives.ciphers import Cipher, algorithms, modes


def validate_key(key: str):
    if not key.isascii():
        raise ValueError("Key must be ASCII-only.")


def decrypt_value_with_prepended_iv_aes_gcm(encrypted_value_base64, key):
    validate_key(key)
    bytes_key = key.encode("utf-8")

    # Decode the base64 encoded encrypted value
    encrypted_value_bytes = base64.b64decode(encrypted_value_base64)

    # Extract the IV (first 12 bytes for AES-GCM)
    iv = encrypted_value_bytes[:12]

    # Extract the ciphertext (excluding the IV and the last 16 bytes which are the authentication tag)
    ciphertext = encrypted_value_bytes[12:-16]

    # Extract the authentication tag (last 16 bytes)
    auth_tag = encrypted_value_bytes[-16:]

    # Create a Cipher object using the extracted IV and the authentication tag
    cipher = Cipher(algorithms.AES(bytes_key), modes.GCM(iv, auth_tag), backend=default_backend())
    decryptor = cipher.decryptor()

    # Decrypt the value
    decrypted_value = decryptor.update(ciphertext) + decryptor.finalize()

    return decrypted_value.decode("utf-8")
