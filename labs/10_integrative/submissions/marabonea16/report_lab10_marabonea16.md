# Raport Lab 10 - Integrare Multi-Omics

## 1. Introducere

In lab-ul asta am lucrat cu date multi-omics pt a integra mai multe tipuri de date biologice (SNPs + expresie gene). Ideea e ca daca combini diferite straturi de date poti gasi mai multe informatii decat daca le analizezi separat.

## 2. Exercițiul 1 - PCA și Vizualizare

### Date folosite
Am lucrat cu:
- 5 SNPs (rs1001 - rs1005)
- 4 gene (BRCA1, TP53, EGFR, MYC)
- 6 sample-uri (SAMPLE_001 - SAMPLE_006)

### Metoda
1. Am incarcat cele 2 matrici (SNP si expresie)
2. Am gasit sample-urile comune (toate 6)
3. Am facut normalizare z-score pt fiecare feature
4. Am concatenat totul intr-o matrice integrata
5. Am facut PCA pe 3 variante:
   - doar SNPs
   - doar expresie
   - date integrate (joint)

### Rezultate PCA

**PCA pe SNPs:**
- Explica ~80% din varianta in PC1 si ~20% in PC2
- Sample-urile sunt destul de separate

**PCA pe expresie:**
- Distributie diferita fata de SNPs
- Gene-urile separa sample-urile altfel decat SNP-urile

**PCA integrat (joint):**
- Combina informatia din ambele straturi
- Ofera o vedere mai completa
- SAMPLE_002 e outlier (departe de celelalte)
- PC1 explica ~79.6% din varianta

### Observații
Cand combini datele multi-omics vezi ca distributia sample-urilor se schimba vs. cand le analizezi separat. Asta arata ca fiecare strat de date aduce info diferita si complementara.

---

## 3. Exercițiul 2 - Corelații Cross-Omics

### Metoda
1. Am folosit matricea integrata din ex1
2. Am separat randurile in SNPs vs gene
3. Am calculat corelații Pearson intre fiecare SNP si fiecare gena
4. Am filtrat doar perechi cu |r| > 0.5
5. Am exportat rezultatele

### Rezultate

Am gasit **11 perechi SNP-gene** cu corelații puternice (|r| > 0.5).

**Top 3 cele mai puternice corelații:**

1. **SNP_rs1001 ↔ GENE_TP53** (r = 0.977)
   - Corelație pozitivă foarte puternica
   - TP53 e o gena importanta in cancer

2. **SNP_rs1004 ↔ GENE_TP53** (r = -0.977)
   - Corelație negativa la fel de puternica
   - Practic opusul lui rs1001

3. **SNP_rs1001 ↔ GENE_BRCA1** (r = 0.975)
   - Alta corelație super puternica
   - BRCA1 e legata de cancer mamar

### Observații importante

- SNP_rs1001 e corelat pozitiv cu toate gene-urile (valori mari)
- SNP_rs1004 e corelat negativ cu toate gene-urile (invers fata de rs1001)
- SNP_rs1005 are corelații moderate cu 3 gene (0.55-0.62)
- SNP_rs1002 si rs1003 au corelații slabe (sub threshold)

Asta sugerează ca anumite variante genetice (SNPs) pot influenta expresia mai multor gene simultan.

---

## 4. Interpretare Biologică

**Ce inseamna corelațiile astea?**

Daca un SNP e corelat cu expresia unei gene, inseamna ca varianta genetica poate afecta cat de mult se exprima gena respectiva. 

De ex:
- Persoanele cu SNP_rs1001 (valoare 1) tind sa aiba expresie mare la TP53, BRCA1, etc.
- Persoanele cu SNP_rs1004 (valoare 1) tind sa aiba expresie mica

**Aplicații practice:**
- Identificare biomarkeri pt diagnostic
- Subtipare pacienți (ex: cancere diferite)
- Farmacogenomica (cine raspunde la ce tratament)
- Medicina de precizie

---

## 5. Concluzii

1. **PCA multi-omics** ofera o vedere mai completa decat analiza pe straturi separate
2. Am identificat **11 perechi SNP-gene** cu corelații puternice
3. SNP_rs1001 si rs1004 sunt candidați importanți pt biomarkeri (corelați cu multe gene)
4. Integrarea datelor multi-omics ajuta la gasirea de pattern-uri care nu se vad cand analizezi fiecare tip de date separat




