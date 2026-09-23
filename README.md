# La crisi della replicazione della scienza: un nuovo approccio statistico bayesiano

Repository della tesi, con codice R per le analisi e per la generazione dei grafici.

## Contenuto della tesi

Il lavoro affronta la *crisi della replicazione* in ambito scientifico, confrontando i metodi classici (frequentisti) con approcci bayesiani per valutare il "successo" di una replicazione. In particolare:

- **Capitolo 1** — Crisi della Replicabilità: cause e possibili soluzioni 
- **Capitolo 2** — Metodi frequentisti (Statistical Significance Criterion, Meta-Analysis SSC, Prediction Interval Criterion, Small Telescope Approach) e metodi bayesiani (Default Bayes Factor, Replication Bayes Factor, Sceptical Bayes Factor, Sceptical Mixture Bayes Factor, Sceptical p-value) applicati ai dati del progetto `RProjects` (pacchetto `ReplicationSuccess`).
- **Capitolo 3** — Calibrazione del parametro γ (grado di "shrinkage" tra studio originale e replica) con tre approcci: calibrazione KL, Empirical Bayes, Full Bayes (modello Stan con diverse prior Beta su γ), incluse diagnostiche MCMC e soglie di classificazione.
- **Capitolo 4** — Simulazioni prior/posterior predittive per un sottoinsieme di studi selezionati, a confronto tra prior informativa (Beta(5,2)) e non informativa (U(0,1)).

## Struttura del repository

```
.
├── README.md
├── tesi/
│   └── tesi_nomecognome.pdf         
├── R/
│   ├── Chapter_2_3_4.R              
│   ├── Plots_Thesis.R               
│   └── Mod_Gamma.stan              
└── .gitignore
```

> **Nota:** lo script `Chapter_2_3_4.R` richiede il file sorgente `Mod_Gamma.stan` nella stessa cartella (viene compilato da `cmdstan_model()`).

## Requisiti

- R (versione ≥ 4.x consigliata)
