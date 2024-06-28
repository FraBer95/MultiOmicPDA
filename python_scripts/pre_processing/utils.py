import pandas as pd


def join_data(rad_path, clin_path, IV_stages):
    df_rad = pd.read_csv(rad_path)
    df_clinical = pd.read_csv(clin_path)

    df_labels = df_clinical[['patient_id', 'df_clinical.tumor_grade_reclass', 'KRAS', 'KRAS_TP53', 'df_clinical.ajcc_pathologic_stage_reclass_2', 'df_clinical.residual_disease_reclass', 'df_clinical.OS_status']]
    df_labels = df_labels.rename(columns={"df_clinical.tumor_grade_reclass" : "tumor_grade_reclass", "df_clinical.ajcc_pathologic_stage_reclass_2" : "pathologic_stage_reclass_2" ,
                                          "df_clinical.residual_disease_reclass" : "residual_disease_reclass", 'df_clinical.OS_status': "status"})

    df_rad = df_rad.rename(columns={"Barcode": "patient_id", "df_clinical.ajcc_pathologic_stage_reclass" : "tumor_grade" })
    # print(df_rad['patient_id'], df_rad['patient_id'].dtypes)
    # print(df_labels['patient_id'], df_labels['patient_id'].dtypes)

    print("Dimensions of clinicals data: ", df_labels.shape)
    print("Dimensions of radiomics data: ", df_rad.shape)
    print("Number of radiomics patients: ", df_rad['patient_id'].value_counts())
    df_merged = pd.merge(df_rad, df_labels, on="patient_id")
    df_total=df_merged.dropna()
    print("Dimensions after inner join", df_total.shape)
    print("Status value count:", df_total['KRAS_TP53'].value_counts())

    conversions = {}
    cat_list = ["pathologic_stage_reclass_2", "residual_disease_reclass", "tumor_grade_reclass"]

    # Iterazione attraverso le colonne categoriche
    for feature in cat_list:
        # Applica la conversione e ottieni il mapping
        encoded, mapping = pd.factorize(df_total[feature])

        # Salva il mapping nella variabile delle conversioni
        conversions[feature] = {'mapping': mapping}

        # Aggiungi la colonna numerica al DataFrame
        df_total[feature] = encoded

    conv = pd.DataFrame(conversions)
    conv.to_csv('../data/conversions.csv', index=False)
    df_no_dup = df_total.drop_duplicates(subset=['patient_id'])
    print("Values count without dup:", df_no_dup['KRAS_TP53'].value_counts())



    print("Dimensions after deleting non relevant features", df_no_dup.shape)


    #new_df.to_csv('../data/matched_dataset.csv', index=False, header=True)
    if not IV_stages: df_no_dup = df_no_dup.drop(df_no_dup[(df_no_dup['tumor_grade'] == 'Stage IV')].index)


    return df_no_dup, df_total


def add_duplicates(df, ids):
    mask = df['patient_id'].isin(ids) #mask based on ids
    new_df = df.loc[mask]   #rows selection based on mask
    new_ids = new_df['patient_id'] #saving ids

    new_df = new_df.iloc[:, 41:] #remove unuseful features,
    concat_df = pd.concat([new_ids, new_df], axis=1)
    return concat_df

def merge_datasets(rad_df, clin_df, gen_df, IVStages):
    if not IVStages:
        ids = ["C3N-03426", "C3N-03430", "C3N-3000", "C3N-01716"]
        clin_df = clin_df[~clin_df['patient_id'].isin(ids)]
        gen_df = gen_df[~gen_df['patient_id'].isin(ids)]
    #values_rad = rad_df['status'].value_counts()
    values_clin = clin_df['status'].value_counts()
    values_mut = gen_df['status'].value_counts()
    clin_df.drop(columns=['status'], inplace=True)
    #gen_df.drop(columns=['status'], inplace=True)

    rad_clin = pd.merge(rad_df, clin_df, on="patient_id")
    clin_mut = pd.merge(clin_df, gen_df, on="patient_id")
    full_data = pd.merge(rad_clin, gen_df, on="patient_id")


   # values = full_data['status'].value_counts()
    # print(f"Status of merged data {values_rad}")
    # print(f"Status of merged data {values_clin}")
    # print(f"Status of merged data {values_mut}")
    # print(f"Status of merged data {values}")
    #df = full_data.drop(columns=["patient_id"])
    return clin_mut




