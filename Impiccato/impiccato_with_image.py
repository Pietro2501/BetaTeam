import os
import random
import string
import tkinter as tk

def pesca_parola(file_words: str):
    """
    Dando in input il percorso ad un file contenente un db di parole, la funzione ne restituisce una a caso e la sua lunghezza
    :type file_word: str
    :param file_word: percorso al file contenente le parole da cui estrarne una casualmente
    :return: tuple
    """
    with open(file_words, "r") as words_db:
        words = words_db.readlines()
        words_clean = [word.strip() for word in words]
        word = random.choice(words_clean)
        words_clean.remove(word)
        return word, len(word), ["_"] * len(word)

def game_control(game: int, game_contest: str = 'choice'):
    """
    La funzione gestisce l'input dell'utente e in base a quello sceglie se procedere nello script o far apparire un messaggio d'errore

    :type game: int
    :type game_contest: string
    :param game: valore intero che, in base al valore assunto, indirizza correttamente il comportamento della funzione
    :param game_contest: stringa che permette sempre di indirizzare il comportamento della funzione, ma a monte in base
                         al contesto in cui la funzione viene chiamata
    :return: int
    """
    if game_contest == 'choice':
        while game != 0 and game != 1:
            try:
                game = int(input("Se vuoi giocare all'IMPICCATO scrivi 1 nell'input box, altrimenti per uscire scrivi 0\n"))
                if game != 0 and game != 1:
                    print("Inserisci solo 0 o 1 come valori numerici")
            except ValueError:
                print("Hai inserito un carattere non numerico")
    elif game_contest == 'mode':
        while game != 0 and game != 1:
            try:
                game = int(input("Se vuoi proporre una lettera contenuta nella parola scrivi 0, se vuoi provare ad indovinare l'intera parola scrivi 1\n"))
                if game != 0 and game != 1:
                    print("Inserisci solo 0 o 1 come valori numerici")
            except ValueError:
                print("Hai inserito un carattere non numerico")
    return game

def load_images():

    images = []

    for i in range(7):
        filename = f"hangman{i}.png"
        try:
            image = tk.PhotoImage(file=filename)
            images.append(image)
        except Exception as e:
            print(f"Errore nel caricamento dell'immagine {filename}: {e}")
            images.append(None)
    return images

def get_image_index(vita):

    image_index = 6 - vita
    if image_index >= len(images):
        image_index = len(images) - 1
    return image_index

def guess_letter(word_list: list, wordwith_: list, game_status: int, vita: int, lista_lettere: list):
    word = "".join(word_list).upper()
    lettera = input('Inserisci una lettera dell\'alfabeto\n')
    lettera = lettera.upper()
    alfabeto = string.ascii_uppercase  # mi restituisce le lettere maiuscole dell'alfabeto
    accenti = "àèìòùéáó"
    # print(alfabeto)
    if alfabeto.find(
            lettera) == -1 and accenti.find(lettera) == -1:  # il .find mi restituisce -1 se la sottostringa specificata (lettera) non è presente nella stringa principale (alfabeto)
        print(f'Hai inserito {lettera}, non una lettera, perdi una vita')
        vita -= 1
        image_label.config(image=images[get_image_index(vita)])
        root.update()
        if vita == 0:
            game_status = 0
            print(f'Hai terminato le vite... la parola da indovinare era {word}\n')
        return wordwith_, vita, game_status
    elif lettera in lista_lettere:
        print(f'Hai già inserito la lettera {lettera}, perdi una vita\n')
        vita -= 1
        image_label.config(image=images[get_image_index(vita)])
        root.update()
        if vita == 0:
            game_status = 0

            print(f'Hai terminato le vite... la parola da indovinare era {word}\n')
        return wordwith_, vita, game_status
    else:
        lista_lettere.append(lettera)

    print(f'Hai inserito la lettera {lettera}\n')
    is_notpresent = True
    for i in range(len(word)):
        if lettera == word[i]:
            wordwith_[i] = lettera
            is_notpresent = False
    if is_notpresent:
        print(f'La lettera {lettera} che hai scelto non è presente nella parola, perdi una vita\n')
        vita -= 1
        image_label.config(image=images[get_image_index(vita)])
        root.update()
        if vita == 0:
            game_status = 0
            print(f'Hai terminato le vite... la parola da indovinare era {word}\n')
        return wordwith_, vita, game_status

    image_label.config(image=images[get_image_index(vita)])
    root.update()

    if wordwith_.count('_') == 0:
        game_status = 0
        print(f'Bravo, hai vinto! La parola era {word}\n')
    return wordwith_, vita, game_status

def guess_word(word_list: list, game_status: bool, vita: int):
    word = "".join(word_list).upper()
    guess_word = input('Indovina la parola: \n')
    guess_word = guess_word.upper()
    if guess_word == word:
        game_status = 0
        print(f'Bravo, hai vinto! La parola era "{word}"\n')
    else:
        print('Hai sbagliato\n')
        vita -= 1

        image_label.config(image=images[get_image_index(vita)])
        root.update()
        if vita == 0:
            game_status = 0

            print(f'Hai terminato le vite... la parola da indovinare era "{word}"\n')

    return game_status, vita


root = tk.Tk()
root.title("Impiccato")


images = load_images()


image_label = tk.Label(root)
image_label.pack()


image_label.config(image=images[0])
root.update()



game_choice = 2
game_choice = game_control(game_choice)

while game_choice == 1:
    lista_letters = []
    parola, len_parola, list_guess = pesca_parola(r"1000_parole_italiane_comuni.txt")
    list_char_word = [i for i in parola.upper()]
    tentativi = 6
    game_stat = 1
    image_label.config(image=images[0])
    root.update()
    print(f"La parola da indovinare è lunga {len_parola} lettere... Cominciamo!")
    print(' '.join(list_guess))
    while game_stat:
        game_mode = 2
        game_mode = game_control(game_mode, 'mode')
        if game_mode:
            game_stat, tentativi = guess_word(list_char_word, game_stat, tentativi)
            if list_guess.count("_") > 0 and game_stat == 1:
                if tentativi > 0:
                    print(f"{' '.join(list_guess)}\nHai ancora {tentativi} tentativi")
        else:
            list_guess, tentativi, game_stat = guess_letter(list_char_word, list_guess, game_stat, tentativi, lista_letters)
            if list_guess.count("_") > 0:
                if tentativi > 0:
                    print(f"{' '.join(list_guess)}\nHai ancora {tentativi} tentativi")
    game_choice = 2
    game_choice = game_control(game_choice)

root.destroy()