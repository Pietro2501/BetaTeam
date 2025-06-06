
# ------------------------------------------- PIETRO DI SPALDRO -------------------------------------------

import os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator
from sklearn.ensemble import RandomForestClassifier
from sklearn.feature_selection import VarianceThreshold
from sklearn.metrics import (accuracy_score, recall_score, confusion_matrix, matthews_corrcoef, roc_auc_score, roc_curve)
from sklearn.model_selection import StratifiedKFold, train_test_split
from sklearn.preprocessing import Normalizer
from sklearn.svm import SVC

# --- Caricamento e preparazione dati ---
df = pd.read_excel("jm101421d_si_001_mod.xls", sheet_name='Dataset', engine='xlrd')
df.columns = df.columns.str.strip() # rimuovo spazi bianchi alle estremità dei nomi delle colonne
smiles = df['SMILES'].astype(str).tolist() # estraggo colonna degli smiles
y_all = df['Activity'].values # estraggo colonna dei label "attività" (0/1)

# --- Generazione Morgan fingerprint (con rdFingerprintGenerator) ---
morgan_generator = rdFingerprintGenerator.GetMorganGenerator(fpSize=512, includeChirality=True)

def morgan_fp_count(smi):
    """
       Genera un fingerprint Morgan count-based a partire da una molecola rappresentata come SMILES.

       Il fingerprint è una rappresentazione numerica della struttura molecolare,
       dove ogni elemento dell'array indica il numero di occorrenze di specifici
       sottostrutture (basate sul raggio e chirality impostati nel generatore Morgan).

       :param smi: Una stringa SMILES che rappresenta la molecola.
       :return: Un array NumPy contenente il fingerprint Morgan (basato sul conteggio delle occorrenze).
       """
    mol = Chem.MolFromSmiles(smi)
    return morgan_generator.GetCountFingerprintAsNumPy(mol)

# Generazione dei fingerprint per tutte le molecole
fps = np.array([morgan_fp_count(s) for s in smiles])
print(f"Ho convertito {len(fps)} molecole in fingerprint array ")

# --- Preprocessing ---
fps = fps[:, fps.sum(axis=0) > 0] # rimozione colonne che sommano a 0
fps = VarianceThreshold(0.02).fit_transform(fps) # eliminazione feature a bassa varianza
X_all = Normalizer().fit_transform(fps) # normalizzazione euclidea

# --- Split train/test (unico seed = 0) ---
idx_all = np.arange(len(X_all)) # recupero gli indici delle molecole
X_tr, X_te, y_tr, y_te, idx_tr, idx_te = train_test_split(X_all, y_all, idx_all, test_size=0.2, random_state=0, stratify=y_all)

# --- Salvataggio fingerprint di train/test ---
pd.DataFrame(X_tr).to_csv("fingerprints_train_seed0.csv", index=False)
pd.DataFrame(X_te).to_csv("fingerprints_test_seed0.csv", index=False)

# --- Salvataggio fingerprint e indici per altri 9 split (seed da 1 a 9) in sotto-cartelle specifiche---
for seed in range(1, 10):
    # creo sotto-cartella per ogni seed
    seed_folder = f"seed_{seed}"
    os.makedirs(seed_folder, exist_ok=True)

    # split
    X_tr_i, X_te_i, idx_tr_i, idx_te_i = train_test_split(X_all, idx_all, test_size=0.2, random_state=seed, stratify=y_all)

    # salvataggio fingerprint e indici nella sotto-cartella
    pd.DataFrame(X_tr_i).to_csv(os.path.join(seed_folder, "fingerprints_train.csv"), index=False)
    pd.DataFrame(X_te_i).to_csv(os.path.join(seed_folder, "fingerprints_test.csv"), index=False)
    pd.DataFrame({'index': idx_tr_i}).to_csv(os.path.join(seed_folder, "indices_train.csv"), index=False)
    pd.DataFrame({'index': idx_te_i}).to_csv(os.path.join(seed_folder, "indices_test.csv"), index=False)

print("\nHo generato correttamente le matrici fingerprint per il seed 0, che userò successivamente per addestrare i modelli.")
print("\nHo generato correttamente 9 cartelle (seed da 1 a 9), ognuna con 4 file csv:"
      "\n- matrice fingerprint delle molecole in train"
      "\n- matrice fingerprint delle molecole in test"
      "\n- colonna di indici delle molecole in train"
      "\n- colonna di indici delle molecole in test\n")


# ------------------------------------------- SALVATORE CHIARIELLO -------------------------------------------

# --- CV interna 5-fold sul training set ---
skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=0) # inizializza la CV
cv_rf_rows, cv_svc_rows = [], [] # liste per salvare predizioni RF e SVC
cv_stats = {'rf': [], 'svc': []} # dizionario per salvare metriche medie

for i, (train_idx, val_idx) in enumerate(skf.split(X_tr, y_tr)):
    X_train, X_val = X_tr[train_idx], X_tr[val_idx] # dati per il fold corrente
    y_train, y_val = y_tr[train_idx], y_tr[val_idx] # label per il fold corrente
    idx_val = idx_tr[val_idx] # indici assoluti delle molecole di validazione

    rf = RandomForestClassifier(random_state=0).fit(X_train, y_train) # addestramento RF
    rf_scores = rf.predict_proba(X_val)[:, 1] # predizione in probabilità
    bin_pred_rf = (rf_scores >= 0.5).astype(int) # predizione label categorica binaria
    tn, fp, fn, tp = confusion_matrix(y_val, bin_pred_rf).ravel() #matrice di confusione
    # memorizzo le performance in cv di quel i-esimo fold in validation
    cv_stats['rf'].append({
        'fold': i + 1,
        'accuracy': round(accuracy_score(y_val, bin_pred_rf),3),
        'sensitivity':  round(recall_score(y_val, bin_pred_rf), 3),
        'specificity':  round(tn / (tn + fp), 3),
        'mcc':  round(matthews_corrcoef(y_val, bin_pred_rf), 3),
        'auc':  round(roc_auc_score(y_val, rf_scores), 3)
    })
    # aggiorno la lista delle predizioni (probabilistiche e categoriche) prodotte da quel i-esimo fold in validation
    cv_rf_rows.extend([{"fold": i + 1, "score": s, "label_true": l, "label_pred": p, "index": int(idx)} for s, l, p, idx in zip(rf_scores, y_val, bin_pred_rf, idx_val)])

    svc = SVC(kernel='linear', probability=True, random_state=0).fit(X_train, y_train) # addestro SVC lineare
    svc_scores = svc.predict_proba(X_val)[:, 1] # predizione in probabilità
    bin_pred_svc = (svc_scores >= 0.5).astype(int) # predizione label categorica binaria
    tn, fp, fn, tp = confusion_matrix(y_val, bin_pred_svc).ravel() # matrice di confusione
    # memorizzo le performance in cv di quel i-esimo fold in validation
    cv_stats['svc'].append({
        'fold': i + 1,
        'accuracy':  round(accuracy_score(y_val, bin_pred_svc), 3),
        'sensitivity':  round(recall_score(y_val, bin_pred_svc), 3),
        'specificity':  round(tn / (tn + fp), 3),
        'mcc':  round(matthews_corrcoef(y_val, bin_pred_svc), 3),
        'auc':  round(roc_auc_score(y_val, svc_scores), 3)
    })
    # aggiorno la lista delle predizioni (probabilistiche e categoriche) prodotte da quel i-esimo fold in validation
    cv_svc_rows.extend([{"fold": i + 1, "score": s, "label_true": l, "label_pred": p, "index": int(idx)} for s, l, p, idx in zip(svc_scores, y_val, bin_pred_svc, idx_val)])

# Salvataggio delle predizioni e metriche in CV per i due modelli
pd.DataFrame(cv_rf_rows).to_csv("cv_rf_all_folds.csv", index=False)
pd.DataFrame(cv_svc_rows).to_csv("cv_svc_all_folds.csv", index=False)
pd.DataFrame(cv_stats['rf']).to_csv("cv_rf_metrics_per_fold.csv", index=False)
pd.DataFrame(cv_stats['svc']).to_csv("cv_svc_metrics_per_fold.csv", index=False)


# ------------------------------------------- CARMELO SMALDINO -------------------------------------------

# --- Test set ---
test_stats = {'rf': {}, 'svc': {}}
rf_full = RandomForestClassifier(random_state=0).fit(X_tr, y_tr) # addestro modello RF su tutto il training set iniziale (80%)
svc_full = SVC(kernel='linear', probability=True, random_state=0).fit(X_tr, y_tr) # addestro modello SVC su tutto il training set iniziale (80%)

test_scr_rf = rf_full.predict_proba(X_te)[:, 1] # score probabilistici RF
test_scr_svc = svc_full.predict_proba(X_te)[:, 1] # score probabilistici SVC

# Dominio di applicabilità
fp_min, fp_max = X_tr.min(axis=0), X_tr.max(axis=0) # definisco due vettori min e max per ogni colonna nel training set
mask_in_domain = np.all((X_te >= fp_min) & (X_te <= fp_max), axis=1) # controllo, per ogni molecola in test, l'appartenenza di ogni feature al range min-max definito in training.
                                                                     # Genero un array booleano in output lungo quanto il numero di molecole in test (true per le molecole che
                                                                     # rispettano l'appartenenza ai range min-max di ogni feature, false per tutte le altre)

# Stampa numero molecole in-domain / out-of-domain
print(f"Test set seed 0: {mask_in_domain.sum()} molecole in-domain, {len(mask_in_domain) - mask_in_domain.sum()} out-of-domain\n")

# definisco la classificazione binaria sulla base della probabilità di predizione
bin_rf = (test_scr_rf >= 0.5).astype(int)
bin_svc = (test_scr_svc >= 0.5).astype(int)
# salvataggio dei risultati in test con indicazione sul dominio
pd.DataFrame({'score': test_scr_rf, 'label_true': y_te, 'label_pred': bin_rf, 'index': idx_te, 'in_domain': mask_in_domain}).to_csv("test_rf.csv", index=False)
pd.DataFrame({'score': test_scr_svc, 'label_true': y_te, 'label_pred': bin_svc, 'index': idx_te, 'in_domain': mask_in_domain}).to_csv("test_svc.csv", index=False)

# Calcolo metriche in test e salvataggio in file di output
for model_name, scores, bin_pred in zip(["rf", "svc"], [test_scr_rf, test_scr_svc], [bin_rf, bin_svc]):
    tn, fp, fn, tp = confusion_matrix(y_te, bin_pred).ravel()
    test_stats[model_name] = {
        'accuracy':  round(accuracy_score(y_te, bin_pred), 3),
        'sensitivity':  round(recall_score(y_te, bin_pred), 3),
        'specificity':  round(tn / (tn + fp), 3),
        'mcc':  round(matthews_corrcoef(y_te, bin_pred), 3),
        'auc':  round(roc_auc_score(y_te, scores), 3)
    }
pd.DataFrame([test_stats['rf']]).to_csv("test_rf_metrics.csv", index=False)
pd.DataFrame([test_stats['svc']]).to_csv("test_svc_metrics.csv", index=False)

def metrics(df):
    """
    Calcola diverse metriche di classificazione a partire da un DataFrame
    contenente le etichette reali, le predizioni binarie e score di predizione del modello.

    :param df: DataFrame Pandas contenente almeno tre colonne:
    :return: Series Pandas con le metriche calcolate, ciascuna arrotondata a 3 cifre decimali.
    """
    tn, fp, fn, tp = confusion_matrix(df['label'], df['pred']).ravel()
    return pd.Series({
        'accuracy': round(accuracy_score(df['label'], df['pred']), 3),
        'sensitivity': round(recall_score(df['label'], df['pred']), 3),
        'specificity': round(tn / (tn + fp), 3),
        'mcc': round(matthews_corrcoef(df['label'], df['pred']), 3),
        'auc': round(roc_auc_score(df['label'], df['score']), 3)
    })

# Metriche separate per in-domain e out-domain in df separati
for model_name, scores, bin_pred in zip(["rf", "svc"], [test_scr_rf, test_scr_svc], [bin_rf, bin_svc]):
    y_pred_df = pd.DataFrame({'score': scores, 'label': y_te, 'pred': bin_pred, 'in_domain': mask_in_domain, 'index': idx_te})
    in_dom = y_pred_df[y_pred_df['in_domain'] == True]
    out_dom = y_pred_df[y_pred_df['in_domain'] == False]

    metrics_in = metrics(in_dom)
    metrics_out = metrics(out_dom)
    metrics_in["group"] = "in_domain"
    metrics_out["group"] = "out_domain"
    combined_metrics = pd.concat([metrics_in, metrics_out], axis=1).T.set_index("group")
    combined_metrics.to_csv(f"{model_name}_domain_performance.csv")

# --- Plot ROC test set completo ---
fpr_rf, tpr_rf, _ = roc_curve(y_te, test_scr_rf)
fpr_svc, tpr_svc, _ = roc_curve(y_te, test_scr_svc)

plt.figure()
plt.plot(fpr_rf, tpr_rf, label=f"RF (AUC={test_stats['rf']['auc']:.3f})")
plt.plot(fpr_svc, tpr_svc, label=f"SVC (AUC={test_stats['svc']['auc']:.3f})")
plt.plot([0, 1], [0, 1], 'k--', label='Random')
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('ROC – Test Set')
plt.legend(loc='lower right')
plt.savefig('roc_test_set.png')

# --- Plot ROC: RF+SVC test set in-domain e out-domain ---
fig, axs = plt.subplots(1, 2, figsize=(14, 6))

# In-domain
auc_rf_in = roc_auc_score(y_te[mask_in_domain], test_scr_rf[mask_in_domain])
auc_svc_in = roc_auc_score(y_te[mask_in_domain], test_scr_svc[mask_in_domain])
fpr_rf_in, tpr_rf_in, _ = roc_curve(y_te[mask_in_domain], test_scr_rf[mask_in_domain])
fpr_svc_in, tpr_svc_in, _ = roc_curve(y_te[mask_in_domain], test_scr_svc[mask_in_domain])
axs[0].plot(fpr_rf_in, tpr_rf_in, label=f"RF (AUC={auc_rf_in:.3f})")
axs[0].plot(fpr_svc_in, tpr_svc_in, label=f"SVC (AUC={auc_svc_in:.3f})")
axs[0].plot([0, 1], [0, 1], 'k--')
axs[0].set_title("In-domain")
axs[0].set_xlabel("False Positive Rate")
axs[0].set_ylabel("True Positive Rate")
axs[0].legend(loc='lower right')

# Out-domain
auc_rf_out = roc_auc_score(y_te[~mask_in_domain], test_scr_rf[~mask_in_domain])
auc_svc_out = roc_auc_score(y_te[~mask_in_domain], test_scr_svc[~mask_in_domain])
fpr_rf_out, tpr_rf_out, _ = roc_curve(y_te[~mask_in_domain], test_scr_rf[~mask_in_domain])
fpr_svc_out, tpr_svc_out, _ = roc_curve(y_te[~mask_in_domain], test_scr_svc[~mask_in_domain])
axs[1].plot(fpr_rf_out, tpr_rf_out, label=f"RF (AUC={auc_rf_out:.3f})")
axs[1].plot(fpr_svc_out, tpr_svc_out, label=f"SVC (AUC={auc_svc_out:.3f})")
axs[1].plot([0, 1], [0, 1], 'k--')
axs[1].set_title("Out-domain")
axs[1].set_xlabel("False Positive Rate")
axs[1].set_ylabel("True Positive Rate")
axs[1].legend(loc='lower right')

plt.tight_layout()
plt.savefig("roc_in_out_domain.png")

print("Completato: \n- Creazione fingerprint\n- Model CV\n- Model Test\n- Model performance\n- Dominio di applicabilità\n- Curve ROC plot")
