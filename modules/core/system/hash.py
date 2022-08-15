from urllib import request
from cryptography.fernet import Fernet
from modules.core.variables.string_man import web_check


class hash:
    def __init__(self, object, id):
        self.object = object
        self.id = id

    def __eq__(self, other):
        return self.object == other.object and self.id == other.id

    def __hash__(self):
        return hash((self.object, self.id))


def generate_key():
    key = Fernet.generate_key()
    with open("secret.key", "wb") as key_file:
        key_file.write(key)


def load_key():
    return open("secret.key", "rb").read()


def encrypt_message(message):
    key = load_key()
    encoded_message = message.encode()
    f = Fernet(key)
    encrypted_message = f.encrypt(encoded_message)

    return encrypted_message


def site_status(
        url,
):
    url = web_check(url)

    if request.urlopen(url).getcode() == 200:
        return True
    else:
        return False
