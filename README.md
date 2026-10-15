# Dengue-patient-classifiers
This repository contains the codes used to analyze whether plasma peptide quantities can be used to classify dengue patients into various disease states and is linked to the medRXIV preprint

Shamorkina, T. M., Nteak, S. K., Lay, S., Kallor, A. A., Ly, S., Duong,V., ... & Snijder, J. (2026). Plasma proteomics identifies early markers of endothelial and inflammatory activation associated with dengue disease severity in children. Medrxiv.

The three classes of disease states are:
```
a) Acute vs Recovery
b) Hospitalized vs Subclinical
c) Dengue Fever (DF) vs Dengue Hemorrhagic Fever (DHF)
```
Plasma peptide LFQ values per patient are normalized, batch-corrected then input to the classifiers, which learns a decision boundary between each class based on these LFQ values.

There are two types of classifiers used:

a) Random forest: This uses a random forest classifier to make decisions about the patient's disease state based on
their plasma LFQ values. Also uses Mean GINI coefficient to determine which genes contribute most significantly to 
the classification task and SHAP calculations to inform what the model is learning.

b) Logistic regression with elastic net: This uses a logistic regression classifier with an elastic net regularization to make decisions about the patient's disease state based on their plasma LFQ values. Also uses Mean GINI coefficient to determine which genes contribute most significantly to the classification task and SHAP calculations to inform what the model is learning.

