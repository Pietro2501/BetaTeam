import os
import random
import string #si tratta di un modulo standard che fornisce strumenti utili per la manipolazione e la gestione delle stringhe


def pesca_parola(file_words: str):
    """
    Dando in input il percorso ad un file contenente una lista di parole, la funzione ne restituisce una a caso e la sua lunghezza

    :type file_word: str
    :param file_word: percorso al file contenente le parole da cui estrarne una casualmente
    :return: tuple
    """
    with open(file_words, "r",encoding="utf-8") as words_db:
        words = words_db.readlines()
        words_clean = [word.strip() for word in words]
        word = random.choice(words_clean)
        words_clean.remove(word)
        return word, len(word), ["_"] * len(word)

def mostraOmino(tentativi):
    """
    Funzione che permette di stampare a video l'omino in base alle vite perse
    
    :type tentativi: int
    :param tentativi: numero di tentativi rimanenti
    
    :return: int
    """
    passi = [
        """
        --------
        |      |
        |      
        |    
        |      
        |     
        -
        """,
        """
        --------
        |      |
        |      O
        |    
        |      
        |     
        -
        """,
        """
        --------
        |      |
        |      O
        |      |
        |      |
        |     
        -
        """,
        """
        --------
        |      |
        |      O
        |     \\|
        |      |
        |     
        -
        """,
        """
        --------
        |      |
        |      O
        |     \\|/
        |      |
        |     
        -
        """,
        """
        --------
        |      |
        |      O
        |     \\|/
        |      |
        |     / 
        -
        """,
        """
        --------
        |      |
        |      O
        |     \\|/
        |      |
        |     / \\
        -
        """
    ]
    return passi[6 - tentativi]

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


def guess_letter(word_list: list, wordwith_: list, game_status: int, vita: int, lista_lettere: list):
    """
    La funzione permette di prendere in input una parola casuale e una lista di '_' lunga quanto la parola. 
    Successivamente gestisce i tentativi di inserire una lettera per indovinare la parola fin quando o si vince o si perde.
    
    :type word_list: list
    :type wordwith_: list
    :type game_status: int
    :type vita: int
    :type lista_lettere: list
    :param word_list: parola selezionata
    :param wordwith_: lista di '_' lunga quanto la parola selezionata
    :param game_status: intero che varia tra 0 (uscire dal gioco) e 1 (continua a giocare) in base a come avanza il gioco
    :param vita: contatore delle vite rimanenti
    :param lista_lettere: lista delle lettere già utilizzate
    
    :return: tuple 
    """
    word = "".join(word_list).upper()
    lettera = input('Inserisci una lettera dell\'alfabeto\n')
    lettera = lettera.upper()
    alfabeto = string.ascii_uppercase  # mi restituisce le lettere maiuscole dell'alfabeto
    accenti = 'ÀÈÌÒÙÉÁÓ'
    # print(alfabeto)
    if alfabeto.find(lettera) == -1 and accenti.find(lettera) == -1:  # il .find mi restituisce -1 se la sottostringa specificata (lettera) non è presente nella stringa principale (alfabeto)
        print(f'Hai inserito {lettera}, non una lettera, perdi una vita')
        vita -= 1
        print(mostraOmino(vita))
        if vita == 0:
            game_status = 0
            print(f'Hai terminato le vite... la parola da indovinare era {word}\n')
        return wordwith_, vita, game_status
    elif lettera in lista_lettere:
        print(f'Hai già inserito la lettera {lettera}, perdi una vita\n')
        vita -= 1
        print(mostraOmino(vita))
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
        print(mostraOmino(vita))
        if vita == 0:
            game_status = 0
            print(f'Hai terminato le vite... la parola da indovinare era {word}\n')
        return wordwith_, vita, game_status

    print(mostraOmino(vita))
    if wordwith_.count('_') == 0:
        game_status = 0
        print(f'Bravo, hai vinto! La parola era {word}\n')
    return wordwith_, vita, game_status

def guess_word(word_list: list, game_status: bool, vita: int):
    """
    La funzione prende in input la parola e permette all'utente di indovinare la parola.
    Se l'utente indovina il gioco finisce, mentre se sbaglia perde una vita e aggiorna l'omino.
    
    :type word_list: list
    :type game_status: bool
    :type vita: int
    :param word_list: parola selezionata
    :param game_status: intero che varia tra 0 (uscire dal gioco) e 1 (continua a giocare) in base a come avanza il gioco
    :param vita: contatore delle vite rimanenti
    
    :return: tuple
    """
    word = "".join(word_list).upper()
    guess_word = input('Indovina la parola: \n')
    guess_word = guess_word.upper()

    if guess_word == word:
        game_status = 0
        print(f'Bravo, hai vinto! La parola era "{word}"\n')
    else:
        print('Hai sbagliato, perdi una vita\n')
        vita -= 1


    print(mostraOmino(vita))

    if vita == 0 and game_status != 0:
        game_status = 0
        print(f'Hai terminato le vite... la parola da indovinare era "{word}"\n')

    return game_status, vita

game_choice = 2
game_choice = game_control(game_choice)

while game_choice == 1:
    lista_letters = []
    parola, len_parola, list_guess = pesca_parola(r"Impiccato\1000_parole_italiane_comuni.txt")
    list_char_word = [i for i in parola.upper()]
    tentativi = 6
    game_stat = 1
    print(f"La parola da indovinare è lunga {len_parola} lettere... Cominciamo!")
    print(mostraOmino(tentativi))
    print(' '.join(list_guess))
    while game_stat:
        game_mode = 2
        game_mode = game_control(game_mode, 'mode')
        if game_mode:
            game_stat, tentativi = guess_word(list_char_word, game_stat, tentativi)
            if game_stat == 1:
                if tentativi > 0:
                    print(f"{' '.join(list_guess)}\nHai ancora {tentativi} tentativi")
        else:
            list_guess, tentativi, game_stat = guess_letter(list_char_word, list_guess, game_stat, tentativi, lista_letters)
            if game_stat == 1:
                if tentativi > 0:
                    print(f"{' '.join(list_guess)}\nHai ancora {tentativi} tentativi")
    game_choice = 2
    game_choice = game_control(game_choice)