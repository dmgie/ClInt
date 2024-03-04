## Functions for vectorized variant analysis
## Use pickle for binary data storage
import pickle

def create_patient_vectors(variants_dict_complete, variants_dict_per_patient):
    
    patient_vectors = {}

    for condition, variants in variants_dict_complete.items():
        patient_vectors[condition] = {}  

        if condition in variants_dict_per_patient:
            for patient, patient_data in variants_dict_per_patient[condition].items():
                patient_vectors[condition][patient] = {}
                for variante in variants:
                    if variante in patient_data:
                        patient_vectors[condition][patient][variante] = patient_data[variante]
                    else:
                        patient_vectors[condition][patient][variante] = "/"

    return patient_vectors                 

def save_patient_vectors(patient_vectors, out_dir):
    with open(f"{out_dir}/vartable_output/patient_vectors.pkl", "wb") as out_file:
        pickle.dump(patient_vectors, out_file)

def save_variant_vectors(variant_vectors, out_dir):
    with open(f"{out_dir}/vartable_output/variant_vectors.pkl", "wb") as out_file:
        pickle.dump(variant_vectors, out_file)
        
def save_runtime(runtimes, out_dir):
    with open(f"{out_dir}/vartable_output/runtimes.pkl", "wb") as out_file:
        pickle.dump(runtimes, out_file)
