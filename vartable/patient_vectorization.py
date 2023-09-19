## Functions for vectorized variant analysis


def create_patient_vectors(variants_dict_complete, variants_dict_per_patient):
    
    patient_vectors = {}
    
    for patient in variants_dict_per_patient:
        patient_vectors[patient] = {}
        for variant_position_complete in variants_dict_complete:
            variant_position = variant_position_complete
            
            if variant_position in variants_dict_per_patient[patient]:
                patient_vectors[patient][variant_position] = variants_dict_per_patient[patient][variant_position]
            else:
                patient_vectors[patient][variant_position] = []
    
    print("\n\n")
    for patient in patient_vectors:
        print(patient, patient_vectors[patient])

    return patient_vectors


def get_position_distance(patient_a, patient_b, position, range):

    distance = 0

    if patient_a[position] == [] or patient_b[position] == []:
        return 1
    else:

        ## Nomianl distance
        if patient_a[position]["variant"] == patient_b[position]["variant"]:
            distance += 0
        else:
            distance += 1
        
        ## Numerical distance
        distance += round((abs(patient_a[position]["counts"] - patient_b[position]["counts"]) / range), 3)

    return distance
        
def get_max_min_counts(data):
    max_counts_per_position = {}
    min_counts_per_position = {}

    for patient, patient_data in data.items():
        for position, values in patient_data.items():
            if values and 'counts' in values:
                counts = values['counts']
                counts_list = counts.split(',')
                for count in counts_list:
                    _, count_value = count.split(':')
                    count_value = int(count_value)
                    if position not in max_counts_per_position or count_value > max_counts_per_position[position]:
                        max_counts_per_position[position] = count_value
                    if position not in min_counts_per_position or count_value < min_counts_per_position[position]:
                        min_counts_per_position[position] = count_value

    return max_counts_per_position, min_counts_per_position