import socket
import webbrowser
import wikipedia
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


def open_web(
        url,
):
    url = sm.web_check(url)
    webbrowser.open(url)
    return


def wiki_search(
        query,
        sound=True,
        text=True,
        sent_num=2,
):
    query = sm.valid_string(query)
    results = wikipedia.summary(query, sentences=sent_num)

    if text:
        print(results)

    return results


def email_strip(
        email,
):
    email = sm.valid_string(email)

    if '@' in email:
        return email[:email.index('@')], email[email.index('@') + 1:]
    else:
        raise ValueError(
            f"Error: Invalid email provided, {email}."
        )


