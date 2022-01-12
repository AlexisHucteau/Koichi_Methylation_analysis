# Koichi_Methylation_analysis
Analysis of DNA methylation through IDHm inhibitor treatment

Git des analyses de methylation de l'ADN en fonction du traitement à l'inhibiteur d'IDHm.

[Back to PhD Project](https://alexishucteau.github.io/PhD_project)

# Koichi Patients analysis

L'idée est de faire 5 groupes :

* Baseline Bon Responder
* Baseline Mauvais Responder
* Post Treatment Bon Responder
* Post Treatment Mauvais Responder
* CD34+

Et de discriminer les mutations IDH1 et IDH2

Comparer Bon et Mauvais responder at Baseline -> DMR_Response

Comparer Baseline et CD34+ -> DMR IDH

Prendre DMR specific de la réponse

Comparer Bon et Mauvais responder at post treatment -> DMR residuel

Comparer Post et Baseline -> DMR traitement

Prendre (DMR residuel - DMR traitement) union avec response specific DMR

Comparer ces response specific DMR avec les IDH specific DMR des WGBS.

# Wilson WGBS analysis

4078 datafiles

| AMPLICON | ATAC-seq | Bisulfite-Seq | ChiP-Seq | OTHER | RNA-Seq | Targeted-Capture | WGS | WXS |
| --------- | --------- | --------- | --------- | --------- | --------- | ---------| ---------| --------- |
| 85 | 39 | 91 | 39 | 33 | 205 | 429 | 2348 | 809 |

Usefull datasets --> Bisulfite-Seq
