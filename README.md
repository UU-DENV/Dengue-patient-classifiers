# Dengue-patient-classifiers
This repository contains the codes used to analyze whether plasma peptide biomarkers can be used to classify dengue patients into severe/subclinical disease states and is linked to the medRXIV preprint

```
Shamorkina, T. M., Nteak, S. K., Lay, S., Kallor, A. A., Ly, S., Duong, V., ... & Snijder, J. (2026). Plasma proteomics identifies early markers of endothelial and inflammatory activation associated with dengue disease severity in children. Medrxiv.
```
There are two types of classifiers used:

a) Random forest: This uses a random forest classifier to make decisions about the patient's disease state based on
their plasma LFQ values. Also uses Mean GINI coefficient to determine which genes contribute most significantly to 
the classification task and SHAP calculations to inform what the model is learning.

b) Logistic regression with elastic net: This uses a logistic regression classifier with an elastic net regularization to make decisions about the patient's disease state based on their plasma LFQ values. Also uses Mean GINI coefficient to determine which genes contribute most significantly to the classification task and SHAP calculations to inform what the model is learning.

