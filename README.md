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

| Pheno A | Pheno B | Number of DMRs | Number of genes | R object |
| ----- | -----| -----| ----- | ----- |
| Bad responder Baseline | Control CD34+ | 131 | 119 | DMR_Bad_Baseline_vs_Control_annotated |
| Good responder Baseline | Control CD34+ | 197 | 127 | DMR_Good_Baseline_vs_Control_annotated |
| Bad responder Baseline | Good responder Baseline | 240 | 321 | DMR_Good_vs_Bad_Baseline_annotated |
| Good Baseline vs Control | Bad Baseline vs Control | 62 | 94 | Specific_Bad_response_annotated |
| Bad responder Post treatment | Control CD34+ | 235 | 185 | DMR_Bad_Post_vs_Control_annotated |
| Good responder Post treatment | Control CD34+ | 109 | 68 | DMR_Good_Post_vs_Control_annotated |
| Good Post vs Control | Bad Post vs Control | 95 | 112 | Specific_Bad_response_post_annotated |


# Wilson WGBS analysis

DMR IDH1 vs other AML subtypes

4388 DMRs

DMR IDH2 vs other AML subtypes

2552 DMRs

# Overlap between Wilson WGBS and Illumina EPIC IDHm inhibitor response

| DMR analysis | number of overlap with IDH1 | number of overlap with IDH2  | parameter changed of overlap | number of overlap with IDH1 | number of overlap with IDH2  |
|---| ---- | ---- |---- | ---- | ----- |
| DMR_Bad_Baseline_vs_Control_annotated | 4 | 2 | max gap = 100 | 5 | 3 |
| DMR_Good_Baseline_vs_Control_annotated | 5 | 2 | max gap = 100 | 5 | 3 |
| DMR_Bad_Posttreatment_vs_Control_annotated |
| Specific_Bad_response_annotated | 0 | 0 | max gap = 100 | 1 | 1 |
|
