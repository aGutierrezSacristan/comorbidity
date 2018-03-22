#!/soft/devel/python-2.7/bin/python
import sys
import datetime

def add_elem_dictionary(dictionary, key, elem, repet = False):
    if key in dictionary:
        aux = dictionary.get(key)
        if repet:
            aux.append(elem)
            dictionary[key] = aux
        else:   
            if not elem in aux:
                aux.append(elem)
                dictionary[key] = aux
    else:
        dictionary[key] = [elem]
    return dictionary

def add_one_dictionary(dictionary, key):

    if key in dictionary:
        aux = dictionary.get(key)
        dictionary[key] = aux + 1
    else:
        dictionary[key] = 1
    return dictionary

def col_num_selection(label_list, header_line):
    col_num_dict =  {}
    i = 0
    while i < len(header_line):
        for lbl in label_list:
            if header_line[i] == lbl:
                col_num_dict[lbl] = i
                label_list.remove(lbl)
                break
        i+=1
    return col_num_dict


def get_dictionary_info(input_path, key_list, value_list, many_values):
    
    column_list = list(set(key_list+value_list))
    
    result_dict = {}
    header = 1
    error_line = {}
    for lin in file(input_path):
        fields = lin.strip().split("\t")
        if header:
            header=0
            col_num_dict = col_num_selection(column_list, fields)
            print col_num_dict
            continue
        """
        if lin in error_line:
            print lin
            #sys.exit()
        error_line[lin] = 1
        """
        key_val_list = []
        for kl in key_list:
            key_val_list.append(fields[col_num_dict[kl]])
        key = "-".join(key_val_list)
        if not len(value_list):
            value = 1
        else:
            val_val_list = []
            for kl in value_list:
                val_val_list.append(fields[col_num_dict[kl]])
            value = "###".join(val_val_list)
        if many_values:
            add_elem_dictionary(result_dict, key, value)
        else:
            result_dict[key] = value
    return result_dict



def run(admis_data_path, codes_path, diag_data_path, patient_data_path, sep_admis_date, sep_birth_date, output_folder_path, first_admin_date, intraCodes = False, aggregatedDis = False):
    print "pat_adm_dict_dates"
    pat_adm_dict_dates= get_dictionary_info(admis_data_path, ["patient_id", "admission_id"], ["admissionStartDate"], False)
    print 2
    pat_gender_dict_birth = get_dictionary_info(patient_data_path, ["patient_id"], ["patient_sex", "patient_dateBirth"], False)
    print 3
    pat_adm_dict_codes = get_dictionary_info(diag_data_path, ["patient_id", "admission_id"], ["diagnoses_code"], True)
    print 4
    
    if aggregatedDis:
        codes_dict = get_dictionary_info(codes_path, ["Code"], ["Agg"], False)
    else:
        codes_dict = get_dictionary_info(codes_path, ["Code"], [], False)
    print 5
    pat_dict_codes = get_dictionary_info(diag_data_path, ["patient_id"], ["diagnoses_code"], True)
    print 6
    pat_in_codes = {}
    for pat_id, codes in pat_dict_codes.iteritems():
        for code in codes:
            if code in codes_dict:
                pat_in_codes[pat_id] = 1
                break
    
    patient_code_dict_year = {}
    patient_code_dict_year_first_date_tomaya = {}
    
    onlyYear = False
    for pat_adm in pat_adm_dict_dates:
        pat_id = "-".join(pat_adm.split("-")[:-1])
        adm_date = pat_adm_dict_dates[pat_adm]
        year_adm_date = adm_date.split(sep_admis_date)[0]
        if len(adm_date.split(sep_admis_date)) > 1:
            adm_format_date = datetime.date(int(year_adm_date), int(adm_date.split(sep_admis_date)[1]), int(adm_date.split(sep_admis_date)[2]))
        else:
            onlyYear = True
            adm_format_date = datetime.date(int(year_adm_date), 1,1)
        diag_codes = pat_adm_dict_codes[pat_adm]
        for code in diag_codes:
            pat_key = pat_id +"-"+code
            if not first_admin_date:
                add_elem_dictionary(patient_code_dict_year, pat_key, adm_format_date)
                
                if pat_key in patient_code_dict_year_first_date_tomaya:
                    aux = patient_code_dict_year_first_date_tomaya[pat_key]
                    if adm_format_date < aux:
                        patient_code_dict_year_first_date_tomaya[pat_key] = adm_format_date
                else:
                    patient_code_dict_year_first_date_tomaya[pat_key] = adm_format_date
                continue
            
            if pat_key in patient_code_dict_year:
                aux = patient_code_dict_year[pat_key]
                if adm_format_date < aux:
                    patient_code_dict_year[pat_key] = adm_format_date
            else:
                patient_code_dict_year[pat_key] = adm_format_date
    
    ofile_all = open(output_folder_path + "all", "w")
    ofile_direc = open(output_folder_path + "direct", "w")
    ofile_final = open(output_folder_path + "final", "w")
    print 7
    lin_error = {}
    
    direct_flag_dict = {}
    for pat_adm in pat_adm_dict_dates:
        pat_id = "-".join(pat_adm.split("-")[:-1])
        #adm_id = pat_adm.split("-")[1]
        #adm_date = pat_adm_dict_dates[pat_adm]
        #year_adm_date = adm_date.split(sep_admis_date)[0]
        diag_codes = pat_adm_dict_codes[pat_adm]
        gender_year = pat_gender_dict_birth[pat_id]
        gender = gender_year.split("###")[0]
        birth_date = gender_year.split("###")[1]
        year_birth = birth_date.split(sep_birth_date)[0]
        
        for code in diag_codes:
            #new_line = []
            #new_line.append(pat_id)
            
            pat_key = pat_id +"-"+code
            admi_final_date_list = []
            if first_admin_date:
                admi_final_date_list.append(patient_code_dict_year[pat_key])
            else:
                admi_final_date_list = patient_code_dict_year[pat_key]
            
            for admi_final_date in admi_final_date_list:
                new_line = []
                new_line.append(pat_id)
                new_line.append(str(admi_final_date.year))
                if intraCodes:
                    YES=False
                    if code in codes_dict:
                        YES = True
                if aggregatedDis:
                    code = codes_dict.get(code, code)
                new_line.append(code)
                new_line.append(gender)
                new_line.append(year_birth)
                new_line.append(str(int(admi_final_date.year)-int(year_birth)))
                
                if "\t".join(new_line) in lin_error:
                    continue
                lin_error["\t".join(new_line)] = 1
                
                ofile_all.write("\t".join(new_line)+ "\n")
                
                    
                if pat_id in pat_in_codes:
                    if intraCodes:
                        if YES:
                            ofile_final.write("\t".join(new_line)+ "\n")
                    else:
                        ofile_final.write("\t".join(new_line)+ "\n")
                    
                    
                    if not first_admin_date:
                        if pat_key in direct_flag_dict:
                            continue
                        direct_flag_dict[pat_key] = 1
                        admi_final_date = patient_code_dict_year_first_date_tomaya[pat_key]
                    
                    
                    if onlyYear:
                        adm_date = admi_final_date.year
                    else:
                        adm_date = sep_admis_date.join([str(admi_final_date.year),str(admi_final_date.month), str(admi_final_date.day)])
                    new_line[1] = str(adm_date)
                    new_line[4] = str(birth_date)
                    ofile_direc.write("\t".join(new_line)+ "\n")
            
    
    ofile_all.close()
    ofile_direc.close() 
    ofile_final.close()

if __name__ == '__main__':
    #python comoRbidity $PATIENT_FOLDER_PATH/ $CODES_FOLDER_PATH/ '-' '-'
    

    PATIENT_FOLDER_PATH = sys.argv[1]
    CODES_FOLDER_PATH = sys.argv[2]
    sep_admis_date = sys.argv[3]
    sep_birth_date = sys.argv[4]
    first_admin_date_param = sys.argv[5]
    intraCodes_param = sys.argv[6]
    aggregatedDis_param = sys.argv[7]

    first_admin_date = False
    if first_admin_date_param == "TRUE":
        first_admin_date = True
        
    intraCodes = False
    if intraCodes_param == "TRUE":
        intraCodes = True
        
        
    aggregatedDis = False
    if aggregatedDis_param == "TRUE":
        aggregatedDis = True

    
 

        
    admis_data_path =PATIENT_FOLDER_PATH + "admissionData.txt"
    codes_path =CODES_FOLDER_PATH +  "indexDiseaseCode.txt"
    diag_data_path =PATIENT_FOLDER_PATH + "diagnosisData.txt"
    patient_data_path = PATIENT_FOLDER_PATH + "patientData.txt"
    
    
    run(admis_data_path, codes_path, diag_data_path, patient_data_path, sep_admis_date, sep_birth_date, PATIENT_FOLDER_PATH, first_admin_date, intraCodes, aggregatedDis)
