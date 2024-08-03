import pandas as pd
import numpy as np
import os
import math
import sys
from sklearn.svm import SVR
from sklearn.model_selection import KFold
from sklearn.metrics import r2_score
from scipy.stats import pearsonr, spearmanr
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import argparse

def get_subjects_info(merge_info_path, mri_file, wgs_file):
    merge_info = pd.read_csv(merge_info_path)  # Info for each participant

    # Subjects who have screening, 6 month, and 12 month
    subject_ids_1y = np.unique([item.split(',')[1].strip('"') for item in open(mri_file).readlines()[1:]])
    wgs_subject_ids = [item.rstrip() for item in open(wgs_file).readlines()]
    subjects_selected = [item for item in wgs_subject_ids if item.upper() in subject_ids_1y]

    # Labels: 0 - Control, 1 - MCI, 2 - Dementia
    labels, labels_bl, labels_m06, labels_m12 = [], [], [], []

    for each in subjects_selected:
        each_upper = each.upper()
        if each_upper == "128_S_1407":
            labels_bl.append(merge_info.loc[(merge_info['PTID'] == each_upper) & (merge_info['VISCODE'] == 'bl')]['DX'].to_list()[0])
            labels_m06.append(merge_info.loc[(merge_info['PTID'] == each_upper) & (merge_info['VISCODE'] == 'm06')]['DX'].to_list()[0])
            labels_m12.append('MCI')
        else:
            labels_bl.append(merge_info.loc[(merge_info['PTID'] == each_upper) & (merge_info['VISCODE'] == 'bl')]['DX'].to_list()[0])
            labels_m06.append(merge_info.loc[(merge_info['PTID'] == each_upper) & (merge_info['VISCODE'] == 'm06')]['DX'].to_list()[0])
            labels_m12.append(merge_info.loc[(merge_info['PTID'] == each_upper) & (merge_info['VISCODE'] == 'm12')]['DX'].to_list()[0])

        subject_info = merge_info.loc[merge_info['PTID'] == each_upper][['VISCODE', 'DX']]
        if 'Dementia' in subject_info['DX'].to_list():
            labels.append(2)
        elif 'MCI' in subject_info['DX'].to_list() and 'Dementia' not in subject_info['DX'].to_list():
            labels.append(1)
        elif 'MCI' not in subject_info['DX'].to_list() and 'Dementia' not in subject_info['DX'].to_list():
            labels.append(0)
        else:
            print('Error')

    return subjects_selected, labels, labels_bl, labels_m06, labels_m12

def get_enformer_feats(subjects_selected, llm_path, chrom, gene_name):
    enformer_feats_p = np.stack([
        np.load(f'{llm_path}/{chrom}_{item}_paternal.npy') for item in subjects_selected
    ])
    enformer_feats_m = np.stack([
        np.load(f'{llm_path}/{chrom}_{item}_maternal.npy') for item in subjects_selected
    ])
    print('Enformer feature:', enformer_feats_p.shape, enformer_feats_m.shape)
    return enformer_feats_p, enformer_feats_m

def integrate_enformer_feats(enformer_feats_p, enformer_feats_m, n_components=7):
    pca_feats_p, pca_feats_m = [], []

    for i in range(enformer_feats_p.shape[1]):
        for j in range(enformer_feats_p.shape[2]):
            pca_feats_p.append(PCA(n_components=n_components).fit_transform(enformer_feats_p[:, i, j, :]))
            pca_feats_m.append(PCA(n_components=n_components).fit_transform(enformer_feats_m[:, i, j, :]))

    pca_feats = np.concatenate(pca_feats_p + pca_feats_m, axis=1)
    return pca_feats

def get_imaging_feats(subjects_selected, img_feat_type):
    df_img_feat = pd.read_csv(f'MRI/parcstats_{img_feat_type}.txt', sep='\t', index_col=0, header=0)
    regions = df_img_feat.columns.to_list()
    img_feat_sc, img_feat_m06, img_feat_m12 = [], [], []

    for each in subjects_selected:
        each_upper = each.upper()
        img_feat_sc.append(df_img_feat.loc[f'{each_upper}_sc'].values)
        img_feat_m06.append(df_img_feat.loc[f'{each_upper}_m06'].values)
        img_feat_m12.append(df_img_feat.loc[f'{each_upper}_m12'].values)

    img_feat_sc = np.stack(img_feat_sc)
    img_feat_m06 = np.stack(img_feat_m06)
    img_feat_m12 = np.stack(img_feat_m12)
    print('Imaging feature:', img_feat_sc.shape, img_feat_m06.shape, img_feat_m12.shape)
    return regions, img_feat_sc, img_feat_m06, img_feat_m12

def get_R_squared(input_feats, img_feat, regions, res_path, apply_norm=True, n_splits=5):
    if apply_norm:
        img_feat = StandardScaler().fit_transform(img_feat)

    with open(res_path, 'w') as f_out:
        for i in range(img_feat.shape[1]):  # Iterate over ROI
            X, y = input_feats, img_feat[:, i]
            kfold_cv = KFold(n_splits=n_splits)
            r2_list, pearsonr_list, spearman_list = [], [], []

            for train_indices, test_indices in kfold_cv.split(X):
                X_train, X_test, y_train, y_test = X[train_indices], X[test_indices], y[train_indices], y[test_indices]
                regressor = SVR()
                regressor.fit(X_train, y_train)
                predictions = regressor.predict(X_test)

                r2_list.append(r2_score(y_test, predictions))
                pearsonr_list.append(pearsonr(y_test, predictions)[0])
                spearman_list.append(spearmanr(y_test, predictions)[0])

            print(f'R2: {np.mean(r2_list):.5f} +/- {np.std(r2_list):.5f}')
            print(f'Pearson r: {np.mean(pearsonr_list):.5f} +/- {np.std(pearsonr_list):.5f}')
            print(f'Spearman r: {np.mean(spearman_list):.5f} +/- {np.std(spearman_list):.5f}')
            f_out.write(f'{regions[i]}\t{np.mean(r2_list)}\t{np.mean(pearsonr_list)}\t{np.mean(spearman_list)}\t{np.std(r2_list)}\t{np.std(pearsonr_list)}\t{np.std(spearman_list)}\n')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate and evaluate Enformer features for a given gene.")
    parser.add_argument('--img_feat_type', type=str, required=True, help='Type of imaging feature.')
    parser.add_argument('--gene_name', type=str, required=True, help='Name of the gene.')
    parser.add_argument('--llm_path', type=str, required=True, help='Path to Enformer feature files.')
    parser.add_argument('--refGene_path', type=str, required=True, help='Path to the refGene hg19 TSS bed file.')
    parser.add_argument('--mri_file', type=str, required=True, help='Path to MRI data file.')
    parser.add_argument('--wgs_file', type=str, required=True, help='Path to WGS subjects file.')
    parser.add_argument('--merge_info_path', type=str, required=True, help='Path to the ADNIMERGE info file.')
    parser.add_argument('--res_path', type=str, required=True, help='Path to save the results.')

    args = parser.parse_args()

    subjects_selected, labels, labels_bl, labels_m06, labels_m12 = get_subjects_info(args.merge_info_path, args.mri_file, args.wgs_file)
    regions, img_feat_sc, img_feat_m06, img_feat_m12 = get_imaging_feats(subjects_selected, img_feat_type=args.img_feat_type)

    gene2loc = {item.split('\t')[4]: (item.split('\t')[0], int(item.split('\t')[1])) for item in open(args.refGene_path).readlines()}
    chrom, center = gene2loc[args.gene_name]

    enformer_feats_p, enformer_feats_m = get_enformer_feats(subjects_selected, llm_path=args.llm_path, chrom=chrom, gene_name=args.gene_name)
    agg_enformer_feats = integrate_enformer_feats(enformer_feats_p, enformer_feats_m)
    get_R_squared(agg_enformer_feats, img_feat_sc, regions, res_path=args.res_path)
