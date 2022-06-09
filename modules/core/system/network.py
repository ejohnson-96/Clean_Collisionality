import socket
import modules.core.variables.string_man as sm

def get_host(

):
    return socket.gethostname()

def get_ip(

):
    host_name = get_host()
    ip_address = socket.gethostbyname(host_name)

    return ip_address

def get_website_ip(
        url,
):
    valid_url = sm.web_check(url)
    return socket.gethostbyname(valid_url)






