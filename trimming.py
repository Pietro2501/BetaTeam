def hard_trimming(quality: list[int], treshold=25):
    quality2 = quality.copy()
    # cicliamo con un for tutti i valori di qualità. Quando troviamo il primo valore < 25,
    # tagliamo tutto ciò che viene dopo e restituiamo una nuova lista fino a quel valore
    for value in range(len(quality2)):
        if quality2[value] < treshold:
            quality2 = quality2[0:value]
            break
    return quality2

def dinamic_trimming(quality: list[int], treshold=25, window=1):
    quality_inv = quality[::-1]
    list_mean = quality_inv[:window]
    if sum(list_mean)/len(list_mean) < treshold:
        for value in range(len(list_mean), len(quality_inv)):
            if (sum(list_mean)+quality_inv[value])/(len(list_mean)+1) < treshold:
                list_mean.append(quality_inv[value])
            else:
                break
        quality_inv_trim = quality_inv[len(list_mean):]
    else: quality_inv_trim = quality_inv
    return quality_inv_trim[::-1]