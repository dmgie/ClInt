# from turtle import shape
from vartable_parser import *
import numpy as np
from sklearn import svm
from sklearn.svm import SVC
from sklearn.model_selection import cross_val_score
from sklearn.metrics import accuracy_score
from sklearn.metrics import classification_report, confusion_matrix
from sklearn.model_selection import cross_val_predict
import matplotlib.pyplot as plt
from sklearn.model_selection import train_test_split
import seaborn as sns

patient_vector, all_variants_list = return_variants_list()
all_variants_list = sorted(all_variants_list)

patient_features_dna = []
patient_features_rna = []
condition_labels = []
counts_feature = []

def feature_importances(coef, names):

    ## Get top 15 highest and lowest SVM weight coefficients

    ## Annotate each feature with either DNA or RNA allele or TPM expression counts 
    fixed_features = ["_DNA_A", "_DNA_C", "_DNA_G", "_DNA_T", "_RNA_A", "_RNA_C", "_RNA_G", "_RNA_T", "_Exp"]
    feature_names = [f"{variant}_{feature}" for variant in names for feature in fixed_features]

    data = list(zip(coef, feature_names))
    sorted_data = sorted(data, key=lambda x: x[0])

    smallest_values = sorted_data[:15]
    largest_values = sorted_data[-15:]

    ## Separate names and values for plotting
    smallest_names, smallest_coeffs = zip(*smallest_values)
    largest_names, largest_coeffs = zip(*largest_values)

    ## SVM weigth coefficient plot
    plt.figure(figsize=(10, 6))
    plt.scatter(smallest_names, smallest_coeffs, label='5 Smallest', color='red')
    plt.scatter(largest_names, largest_coeffs, label='5 Largest', color='blue')
    plt.xlabel('Coefficients')
    plt.ylabel('Feature Names')
    plt.title('Top 10 smallest and largest SVM coefficients')
    plt.savefig('weights.png', dpi=300, bbox_inches='tight')

## Create SVM dataset
#### -> each variant hold 9 features, DNA/RNA alleles (8 characters) and expression (1 character)
#### -> Test COVIDsevere condtion vs. rest (otherwise too many misclassifications)

for condition in patient_vector:
    for patient in patient_vector[condition]:
        patient_variants = patient_vector[condition][patient]
        
        dna_gt_features = []
        rna_gt_features = []
        counts = []

        if condition in ["COVIDsevere"]:
            condition_labels.append(condition)
        else:
            condition_labels.append("REST")
        

        for variant in all_variants_list:
            if variant in patient_variants:
                locus = patient_variants[variant]

                if "SNP" in locus:
                    dna_gt_features.append(locus["SNP"]["dna_nuc"])
                    rna_gt_features.append(locus["SNP"]["rna_nuc"])
                    counts.append(locus["counts"])
                    
                else:
                    dna_gt_features.append([0,0,0,0])
                    rna_gt_features.append([0,0,0,0])
                    counts.append(0)

            else:
                dna_gt_features.append([0,0,0,0])
                rna_gt_features.append([0,0,0,0])
                counts.append(0)
                

        patient_features_dna.append(dna_gt_features)
        patient_features_rna.append(rna_gt_features)
        counts_feature.append(counts)



## Concatenate feature vectors
counts_feature = np.array(counts_feature)[:, :, np.newaxis]
combined_features = np.concatenate((np.array(patient_features_dna), np.array(patient_features_rna)), axis=2)
combined_features_expression = np.concatenate((combined_features, np.array(counts_feature)), axis=2)

## Create SVM features and lables
X = combined_features_expression
X = X.reshape(54, -1)
y = np.array(condition_labels)

## Train/Test split
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size = 0.20)

clf = svm.SVC(C=0.1, gamma=1, kernel='linear')
clf.fit(X, y)

## GridSearch for parapeter fine-tuning

# param_grid = {'C': [0.1,1, 10, 100], 'gamma': [1,0.1,0.01,0.001],'kernel': ['linear']}
# grid = GridSearchCV(SVC(),param_grid,refit=True,verbose=2)
# grid.fit(X_train,y_train)
# print(grid.best_estimator_)

# grid_predictions = grid.predict(X_test)
# print(confusion_matrix(y_test,grid_predictions))
# print(classification_report(y_test,grid_predictions))#Output


## Determine most important features

feature_importances(clf.coef_[0], all_variants_list)
predicted = cross_val_predict(clf, X, y, cv=6)
accuracy = accuracy_score(y, predicted)
print(f'Overall Accuracy: {accuracy:.2f}')

## Calculate confusion matrix
conf_matrix = confusion_matrix(y, predicted)
print("Confusion Matrix:")
print(conf_matrix)

## Plot confusion matrix

plt.figure(figsize=(6, 4))
sns.heatmap(conf_matrix, annot=True, fmt="d", cmap="Blues", linewidths=.5, xticklabels=["severe", "rest"], yticklabels=["severe", "rest"], cbar=False)
plt.title('Confusion Matrix')
plt.xlabel('Predicted')
plt.ylabel('True')
plt.savefig('confusion_matrix.png', dpi=300, bbox_inches='tight')
