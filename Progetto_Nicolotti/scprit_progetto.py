import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import rdFingerprintGenerator
from sklearn.ensemble import RandomForestClassifier
from sklearn.feature_selection import VarianceThreshold
from sklearn.metrics import (
    accuracy_score, recall_score, confusion_matrix,
    matthews_corrcoef, roc_auc_score, roc_curve
)
from sklearn.model_selection import train_test_split, StratifiedKFold
from sklearn.preprocessing import Normalizer
from sklearn.svm import SVR

# --- 0. Disabilito i log di RDKit (warning) ---
RDLogger.DisableLog('rdApp.*')

# --- 1. Caricamento e preparazione dati ---
df = pd.read_excel('jm101421d_si_001.xls', sheet_name='Dataset', engine='xlrd')
df.columns = df.columns.str.strip() # rimuovo spazi bianchi alle estremità dei nomi delle colonne
if 'Name' in df.columns:
    df = df.drop(columns=['Name']) # rimuovo la colonna "name" (inutile)

smiles = df['SMILES'].astype(str).tolist() # estraggo colonna degli smiles
y_all = df['Activity'].values # estraggo colonna dei label "attività" (0/1)

# --- 2. Generazione Morgan fingerprint (con rdFingerprintGenerator) ---

# Generatore Morgan fingerprint (riutilizzabile)
morgan_generator = rdFingerprintGenerator.GetMorganGenerator(fpSize=512,includeChirality=True) # oggetto generatore fp

def morgan_fp(smi):
    """
    Questa funzione genera l'oggetto molecola (mol) dall'oggetto smile; successivamente
    da questa genera l'oggetto fingerprint. Infine, per ottenere il vettore binario
    di nostro interesse, dobbiamo creare un array di "0" con NumPy in cui andremo
    a porre pari a "1" solo gli indici dei bit "attivi" indicati dal'oggetto fp
    :param smi:
    :return: arr (list)
    """
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        return None
    fp = morgan_generator.GetFingerprint(mol)
    arr = np.zeros(512, dtype=int)
    arr[list(fp.GetOnBits())] = 1
    return arr

fps_list = [morgan_fp(s) for s in smiles] # list comprehension degli array binari ottenuti applicando la funzione a tutti gli smiles
print(f"Ho convertito {len(fps_list)} molecole in fingerprint array ")
valid = [i for i,f in enumerate(fps_list) if f is not None] #memorizzo gli indici, da fps_list, delle fp valide
print(f"Ho mantenuto poi {len(valid)} fingerprint array validi")
fps = np.array([fps_list[i] for i in valid]) # genero una nuova lista degli array validi
y_all = y_all[valid] # recupero anche gli elementi dalla colonna delle attività corrispondenti agli indici

# --- 3. Preprocessing ---
fps = fps[:, fps.sum(axis=0) > 0] # rimuovo molecole con bit sempre zero
fps = VarianceThreshold(0.02).fit_transform(fps) # rimuovo molecole con varianza ridotta (<0.02)
X_all = Normalizer().fit_transform(fps) # ogni vettore scalato per avere norma unitaria (evitando che quelli più lunghi abbiano più peso solo per la lunghezza)

# --- 4. CV interna + valutazione test across seeds ---
seeds                = list(range(10))
cv_stats             = {'RF': [], 'SVR': []} # ogni lista conterrà 10 dizionari annidati (per ogni seed), ognuna con le metriche delle performance in CV
test_stats           = {'RF': [], 'SVR': []} # ogni lista conterrà 10 dizionari annidati (per ogni seed), ognuna con le metriche delle performance sul test set
test_scores_per_seed = {'RF': [], 'SVR': []} # memorizzo valori di predizione CONTINUI (non binari), per ogni seed, nel test (utile per AUC e curva ROC)
y_tests_per_seed     = [] # qui avrò liste annidate, per ogni seed, delle label di classificazione binaria predette in test

for seed in seeds:
    # 4.1 split esterno (suddivisone proporzionata tra train e test)
    X_tr, X_te, y_tr, y_te = train_test_split(X_all, y_all,test_size=0.2,random_state=seed,stratify=y_all )

    # 4.2 CV interna 5-fold sul training set
    skf   = StratifiedKFold(n_splits=5, shuffle=True, random_state=seed) # variante cv che mantiene le proporzioni delle classi
    # creo un dizionario con due liste, ognuna con array di 0 inizialmente (utili a memorizzare la predizione binaria dell'attività)
    preds = {'RF': np.zeros(len(y_tr)), 'SVR': np.zeros(len(y_tr))} # ogni indice verrà riempito nel momento in cui quella molecola finirà nel validation
    for train_idx, val_idx in skf.split(X_tr, y_tr): # il metodo split restituisce un generatore che, ad ogni iterazione, produce una coppia di array di indici
        X_train, X_val = X_tr[train_idx], X_tr[val_idx] # creo sub-set correnti di fp binarie, per train e validation, sfruttando gli indici
        y_train, y_val = y_tr[train_idx], y_tr[val_idx] # creo sub-set correnti di label attività "0/1", per train e validation, sfruttando gli indici

        # Random Forest
        rf = RandomForestClassifier(random_state=seed)
        rf.fit(X_train, y_train)
        # uso metodo predict_proba per ottenere le probabilità di appartenenza ad una delle due classificazioni "0/1" (es. [0.3, 0.7])
        # ottengo in output una matrice di dimensioni (n_val_samples, 2)
        preds['RF'][val_idx] = rf.predict_proba(X_val)[:,1] # isolo la colonna di probabilità di appartenenza alla classe "1" e associo
                                                            # i valori agli indici corrispondenti nel dizionario preds[RF]

        # SVR
        svr = SVR(kernel='linear')
        svr.fit(X_train, y_train)
        # essendo una regressione, il metodo predict qui restituisce una valore continuo tra 0 e 1,
        # score che rappresentano la confidenza del modello nel classificare come "1"
        preds['SVR'][val_idx] = svr.predict(X_val) # associo anche qui i valori ai corrispondenti indici in preds[SVR]

    # 4.3 raccolgo metriche CV (validation folds)
    for model in ['RF','SVR']:
        scr = preds[model] # isolo gli score (scr) di predizione in validation per ogni modello
        bin_pred = (scr >= 0.5).astype(int) # creo lista binaria di 0/1, leggendo scr, usando la soglia di riferimento
        tn, fp, fn, tp = confusion_matrix(y_tr, bin_pred).ravel() # confronto le etichette originali con quelle predette in cv
        # associo ora le metriche di predizione in cv al sotto-dizionario (seed-chiave sottinteso) corrispondente
        cv_stats[model].append({
            'accuracy': accuracy_score(y_tr, bin_pred),
            'sensitivity': recall_score(y_tr, bin_pred),
            'specificity': tn/(tn+fp),
            'mcc': matthews_corrcoef(y_tr, bin_pred),
            'auc': roc_auc_score(y_tr, scr)
        })

    # 4.4 valutazione su test set
    rf_full      = RandomForestClassifier(random_state=seed).fit(X_tr, y_tr) # riaddestro sull'intero test set (80% di partenza)
    svr_full     = SVR(kernel='linear').fit(X_tr, y_tr)
    test_scr_rf  = rf_full.predict_proba(X_te)[:,1] # predico sul test set (20% di partenza)
    test_scr_svr = svr_full.predict(X_te)

    # salvo scores di predizione per il seed corrente
    test_scores_per_seed['RF'].append(test_scr_rf)
    test_scores_per_seed['SVR'].append(test_scr_svr)
    # salvo la classificazione predetta per il seed corrente
    y_tests_per_seed.append(y_te.copy())

    # calcolo performance in test per ogni seed
    for model, scr, y_true in [('RF',  test_scr_rf,  y_te), ('SVR', test_scr_svr, y_te)]:
        bin_pred = (scr >= 0.5).astype(int) # creo lista binaria di 0/1, leggendo scr, usando la soglia di riferimento
        tn, fp, fn, tp = confusion_matrix(y_true, bin_pred).ravel() # confronto le etichette originali con quelle predette in test
        # associo ora le metriche di predizione in test al sotto-dizionario (seed-chiave sottinteso) corrispondente
        test_stats[model].append({
            'accuracy': accuracy_score(y_true, bin_pred),
            'sensitivity': recall_score(y_true, bin_pred),
            'specificity': tn/(tn+fp),
            'mcc': matthews_corrcoef(y_true, bin_pred),
            'auc': roc_auc_score(y_true, scr)
        })

    # 4.5 (OPZIONALE!!! e un pò inutile forse...) salva predizioni "seed=0"
    if seed == 0:
        pd.DataFrame({
            'orig_index': np.arange(len(y_tr)),
            'score_RF': preds['RF'],
            'score_SVR': preds['SVR']
        }).to_csv('cv_predictions_seed0.csv', index=False)
        pd.DataFrame({
            'orig_index': np.arange(len(y_te)),
            'score_RF': test_scr_rf,
            'score_SVR': test_scr_svr
        }).to_csv('test_predictions_seed0.csv', index=False)

# --- 5. Esporto summary mean±std in un’unica colonna ---
def make_summary(stats_dict, out_csv):
    """
    Prendo in input il dizionario con le metriche delle performance (per cv o test) e la stringa per il nome del csv
    da restituire in output. Calcolo media e deviazione standard e restituisco un file csv in output.
    :param stats_dict (dict)
    :param out_csv (str)
    :return: dataframe
    """
    rows = []
    for model in ['RF','SVR']:
        dfm = pd.DataFrame(stats_dict[model]) # trasformo i sottodizionari (seed:metriche) associati ad un modello in dataframe
        m, s = dfm.mean(), dfm.std() # calcolo media e dev. standard
        # produco una sotto-lista di caratteristiche per quel modello
        rows.append({
            'Model': model,
            'Accuracy (mean±std)':    f"{m.accuracy:.3f}±{s.accuracy:.3f}",
            'Sensitivity (mean±std)': f"{m.sensitivity:.3f}±{s.sensitivity:.3f}",
            'Specificity (mean±std)': f"{m.specificity:.3f}±{s.specificity:.3f}",
            'MCC (mean±std)':         f"{m.mcc:.3f}±{s.mcc:.3f}",
            'AUC (mean±std)':         f"{m.auc:.3f}±{s.auc:.3f}"
        })
    pd.DataFrame(rows).to_csv(out_csv, index=False, encoding='utf-8-sig') # ogni sotto-lista sarà una riga della tabella prodotta

make_summary(cv_stats,   'cv_performance_summary.csv')
make_summary(test_stats, 'test_performance_summary.csv')

# --- 6. Plot ROC per “best seed” di ciascun modello ---

# calcolo gli AUC score, per i due modelli, iterando sui seed e li salvo in una lista
auc_list_rf  = [roc_auc_score(y_tests_per_seed[i], test_scores_per_seed['RF'][i]) for i in range(len(seeds))]
auc_list_svr = [roc_auc_score(y_tests_per_seed[i], test_scores_per_seed['SVR'][i]) for i in range(len(seeds))]

# isolo il seed associato al miglior AUC score per i due modelli
best_rf_seed  = int(np.argmax(auc_list_rf))
best_svr_seed = int(np.argmax(auc_list_svr))

# calcolo fpr/tpr
fpr_rf,  tpr_rf,  _ = roc_curve(y_tests_per_seed[best_rf_seed],
                                test_scores_per_seed['RF'][best_rf_seed])
fpr_svr, tpr_svr, _ = roc_curve(y_tests_per_seed[best_svr_seed],
                                test_scores_per_seed['SVR'][best_svr_seed])

plt.figure()
plt.plot(fpr_rf,  tpr_rf,  label=f"RF (seed={best_rf_seed}, AUC={auc_list_rf[best_rf_seed]:.2f})")
plt.plot(fpr_svr, tpr_svr, label=f"SVR(seed={best_svr_seed}, AUC={auc_list_svr[best_svr_seed]:.2f})")
plt.plot([0,1],[0,1],'k--', label='Random')
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('ROC – Best Seeds per Model')
plt.legend(loc='lower right')
plt.savefig('roc_best_seeds.png')

print("Generati:")
print("- cv_performance_summary.csv")
print("- test_performance_summary.csv")
print("- cv_predictions_seed0.csv")
print("- test_predictions_seed0.csv")
print("- roc_best_seeds.png")
