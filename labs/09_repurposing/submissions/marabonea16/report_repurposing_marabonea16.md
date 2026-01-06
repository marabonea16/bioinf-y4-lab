
Am folosit un set de date care continea relatii drug-gen (care medicament tinteste care gene) cu urmatoarele caracteristici:

- **Total interactiuni:** 16 perechi drug-gen
- **Medicamente unice:** 8 (Aspirin, Cisplatin, Doxorubicin, Etoposide, Ibuprofen, Methotrexate, Nutlin-3, PARP_inhibitor)
- **Gene unice:** 10 (TP53, BRCA1, MDM2, PTGS1, PTGS2, TOP2A, PARP1, TNF, DHFR, TYMS)
- **Gene asociate bolii:** 4 (TP53, BRCA1, MDM2, PTGS2)

Un lucru interesant e ca fiecare medicament din dataset tinteste exact 2 gene, ceea ce face datele destul de echilibrate pentru analiza.

Am creat un **graf bipartit** G = (V, E) unde:
- **V = V_drug ∪ V_gene** (doua seturi disjuncte de noduri)
- **E** = muchiile care conecteaza medicamentele la genele pe care le tintesc

**Proprietatile retelei:**
- **Noduri:** 18 (8 medicamente + 10 gene)
- **Muchii:** 16 interactiuni drug-gen
- **Tip:** Graf neorientat, bipartit


Pentru a vedea cat de asemanate sunt medicamentele, am folosit **coeficientul Jaccard** pe multimile de gene tinta:

$$J(D_i, D_j) = \frac{|G_i \cap G_j|}{|G_i \cup G_j|}$$

unde $G_i$ si $G_j$ sunt seturile de gene tintite de medicamentele $D_i$ si $D_j$.

**Ce inseamna rezultatele:**
- **J = 1.0:** Medicamentele tintesc exact aceleasi gene (similitudine perfecta)
- **0 < J < 1:** Medicamentele au unele gene in comun (partial asemanate)
- **J = 0:** Niciun gene in comun (mecanisme complet diferite)


Pentru fiecare medicament, am calculat **distanta medie** in retea pana la genele care sunt asociate bolii:

$$d_{avg}(D_i) = \frac{1}{|G_{disease}|} \sum_{g \in G_{disease}} d(D_i, g)$$

unde $d(D_i, g)$ este lungimea celui mai scurt drum in graf entre medicamentul $D_i$ si gena $g$.

**Ce inseamna distantele:**
- **d = 1:** Medicamentul tinteste direct gena bolii
- **d = 2:** Medicamentul tinteste o gena care interactioneaza cu gena bolii (efect indirect)
- **d > 2:** Relatie mai indirecta, potential mai slaba



**Top 5 perechi de medicamente cu cea mai mare similaritate:**

| Medicament 1 | Medicament 2 | Similaritate | Ce am observat |
|--------|--------|------------|----------------|
| Cisplatin | Doxorubicin | 1.00 | Tintesc exact aceleasi gene (TP53, BRCA1) |
| Aspirin | Ibuprofen | 0.33 | Amandoua afecteaza calea PTGS2 |
| Cisplatin | Etoposide | 0.33 | Ambele afecteaza calea p53 |
| Cisplatin | Nutlin-3 | 0.33 | Ambele moduleaza semnalarea p53 |
| Doxorubicin | PARP_inhibitor | 0.33 | Ambele sunt implicate in raspunsul la daunele ADN |

**Observatii importante:**
- **Cisplatin si Doxorubicin** au similaritate perfecta (1.0) - amandoua sunt agenti care deterioreaza ADN si activeaza proteina p53
- Eu observ ca mai multe medicamente se grupeaza in jurul **caii p53** (Cisplatin, Doxorubicin, Nutlin-3, Etoposide) - asta are sens din punct de vedere biologic
- **Medicamentele anti-inflamatorii** (Aspirin, Ibuprofen) formeaza un grup separat care tintesc caile ciclooksigenazei


Am calculat distanta fiecarui medicament pana la genele bolii (TP53, BRCA1, MDM2, PTGS2) si am clasat-o:

| Loc | Medicament | Distanta | Tinteste care gene | Observatii |
|------|------|----------|-------|----------------|
| 1 | Cisplatin | 2.75 | TP53, BRCA1 | Direct in caile cancerului |
| 2 | Doxorubicin | 2.75 | TP53, BRCA1 | Activeaza direct p53 |
| 3 | Nutlin-3 | 2.75 | MDM2, TP53 | Inhibitor MDM2, stabilizeaza p53 |
| 4 | Etoposide | 3.25 | TP53, TOP2A | Activeaza p53 prin deteriorarea ADN |
| 5 | PARP_inhibitor | 3.75 | BRCA1, PARP1 | Pentru cancerele cu mutatii BRCA |
| 6 | Aspirin | 4.75 | PTGS1, PTGS2 | Anti-inflamator, mai departe de genele bolii |
| 7 | Ibuprofen | 4.75 | PTGS2, TNF | Anti-inflamator |
| 8 | Methotrexate | 6.00 | DHFR, TYMS | Metanolism folatelor, mai departe |

**Distributia rezultatelor:**
- **Top tier (d < 3.0):** 3 medicamente care tintesc direct genele bolii
- **Nivel mediu (3.0 ≤ d < 4.0):** 2 medicamente cu conexiuni indirecte
- **Bottom (d ≥ 4.5):** 3 medicamente mai departe de genele bolii



**Candidatii cu cea mai mare prioritate (distanta ≤ 2.75):**

1. **Cisplatin si Doxorubicin** — Medicamente oncologice cunoscute
   - Ambele tintesc direct TP53 si BRCA1
   - Mecanismul: Deteriorare ADN → activare p53 → apoptoza
   - Sunt deja folosite in terapia cancerului
   - Potential de refolosire: Pentru alte tipuri de cancer cu probleme p53/BRCA1

2. **Nutlin-3** — Compus experimental promitator
   - Inhibitor MDM2 (care de obicei degradeaza p53)
   - Mecanismul: Stabilizeaza p53 fara a deteriora ADN (mai putin toxic)
   - Studiat pentru tumori cu p53 normal
   - Potential: Activare p53 mai sigura, fara efecte genotoxice



Analiza m-a ajutat sa identific trei grupe principale de medicamente:

**Grupa 1: Raspunsul la daunele ADN (d = 2.75-3.25)**
- Medicamente: Cisplatin, Doxorubicin, Etoposide, Nutlin-3
- Caia comuna: Activarea p53
- Context terapeutic: Cancer (apoptoza dependenta de p53)

**Grupa 2: Repararea ADN (d = 3.75)**
- Medicament: PARP_inhibitor
- Mecanismul: Sintetica letalitate cu deficienta BRCA1
- Context: Cancere cu mutatii BRCA (ovar, san)

**Grupa 3: Anti-inflamatorii (d = 4.75)**
- Medicamente: Aspirin, Ibuprofen
- Mecanismul: Inhibitie ciclooksigenaza
- Context: Mai departe de genele principale ale bolii, potential rol preventiv



- **Medicamentele cu similaritate mare au de obicei si distante similare** (ex: Cisplatin-Doxorubicin)
- **Grupele biologice apar natural din topologia retelei**
- **Distanta in retea coreleaza cu relevanta terapeutica cunoscuta** -



**Rezultatele principale:**
- **Cisplatin, Doxorubicin si Nutlin-3** sunt cei mai buni candidati pentru context bolii (implicare p53)
- Proximitatea in retea coreleaza cu relevanta terapeutica
- Abordarea ofera o metoda sistematica si bazata pe date pentru refolosire


