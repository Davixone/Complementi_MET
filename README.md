# Maxwell-Zener Viscoelastic Oscillator Dynamics

Questo repository contiene il codice MATLAB sviluppato per l'elaborato di **Complementi di Meccanica e Termodinamica** (A.A. 2025/2026), Università degli Studi dell'Aquila.

Il progetto analizza la dinamica di un oscillatore viscoelastico basato sul modello "Standard Linear Solid" (Maxwell-Zener).

## Contenuto del Repository

Gli script permettono di riprodurre tutte le figure presenti nella tesi:

* **Analisi Stabilità:** `Delta_cubico_3D.m` visualizza il discriminante cubico.
* **Risposte Temporali:** `creep_maxwell_zener.m` e `relaxation_maxwell_zener.m` per le risposte al gradino.
* **Dinamica Forzata:** `forced_maxwell_zener_response.m` calcola la risposta in frequenza $H(i\omega)$ e la risonanza.
* **Spazio delle Fasi:** `phase_space_maxwell_zener.m` traccia le orbite 3D.
* **Validazione:** `comparison_analytic_numerical.m` confronta la soluzione analitica con `ode45`/`ode15s`.

## Requisiti
* MATLAB R2020b o successivo (testato su R2024b).
* Non sono richiesti toolbox aggiuntivi.

## Autore
Davide G. Vigna
