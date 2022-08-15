from PyDictionary import PyDictionary
from spellchecker import SpellChecker
from modules.core.variables.string_man import valid_string
from modules.core.variables.char_man import lower_all_letter, capital_first_letter

dictionary = PyDictionary()


def word_exists(
        word,
        dict_type="en_GB",
):
    return


def spell_check(
        word,
        dict_lang='en',
):
    valid_languages = {'English': 'en', 'Spanish': 'es', 'French': 'fr',
                       'Portuguese': 'pt', 'German': 'de',
                       'Russian': 'ru'}
    dict_lang = lower_all_letter(valid_string(dict_lang))

    valid_language = False
    for language in valid_languages:
        if language == dict_lang:
            dict_lang = valid_languages[language]
            valid_language = True
        if valid_languages[language] == dict_lang:
            valid_language = True

    if not valid_language:
        raise ValueError(
            "Error: Dictionary language provided is not supported, "
            f"provided {dict_lang}, supported languages are: "
            f"{valid_languages.keys()}."
        )

    spell = SpellChecker()
    word = lower_all_letter(valid_string(word))

    if spell.unknown(word) == spell.correction(word):
        return True
    else:
        return False


def spell_review(
        word,
        correction=True,
        suggestion=True,
):
    word = valid_string(word)
    print(spell_check(word))
    if spell_check(word):
        raise ValueError(
            "Error: Argument provided is spelt correctly."
        )
    else:
        spell = SpellChecker()
        if correction and not suggestion:
            return spell.correction(word)
        elif suggestion and not correction:
            return spell.candidates(word)
        elif correction and suggestion:
            return spell.correction(word), spell.candidates(word)
        else:
            raise ValueError(
                "Error: No value returned, one bool must be true."
            )

    return


y = spell_review('science', suggestion=True, correction=True)
print( y)
