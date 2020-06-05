'''   
   <CombinationFinder identifies 4-hit combinations of carcinogenic genes.>
    Copyright (C) <2020>  <CombinationFinderDevs>
    
    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
'''

from statsmodels.stats.proportion import proportion_confint as ci
import numpy as np

def loadGeneSampleMatrixTumor(tumor_file):
    gene = []
    sample = []
    gene2id = {}

    reader = open(tumor_file, "r")
    header = reader.readline()
    tokens = header.split()
    num_genes = int(tokens[0])
    num_samples = int(tokens[1])
    gene_sample_matrix = [[0 for _ in range(num_samples)] for _ in range(num_genes)]

    for line in reader:
        line = line.strip()
        tokens = line.split()
        gene_id, sample_id, mutated, gene = int(tokens[0]), int(tokens[1]), int(tokens[2]), tokens[3]
        if mutated > 0:
            mutated = 1        
        gene_sample_matrix[gene_id][sample_id] = mutated
        gene2id[gene] = gene_id

    reader.close()

    return gene_sample_matrix, gene2id

def loadGeneSampleMatrixNormal(normal_file, gene2id):
    sample_id_map = {}
    current_sample_id = 0
    num_genes = len(gene2id)
 
    with open(normal_file, "r") as reader:
        for line in reader:
            line = line.strip()
            tokens = line.split()
            sample_id = int(tokens[1])
            if sample_id not in sample_id_map:
                sample_id_map[sample_id] = current_sample_id
                current_sample_id += 1

    num_samples = current_sample_id
    gene_sample_matrix = [[0 for _ in range(num_samples)] for _ in range(num_genes)]
    
    with open(normal_file, "r") as reader:
        for line in reader:
            line = line.strip()
            tokens = line.split()
            sample_id = sample_id_map[int(tokens[1])]
            if tokens[0] in gene2id:
                gene_id = gene2id[tokens[0]]
                gene_sample_matrix[gene_id][sample_id] = 1

    return gene_sample_matrix


def get_covered_samples(gene_sample_matrix, gene2id, combination):
    hits = len(combination)
    num_samples = len(gene_sample_matrix[0])
    covered_samples = set(range(num_samples))
    for gene in combination:
        samples = set()
        gid = gene2id[gene]
        for sample in range(num_samples):
            if gene_sample_matrix[gid][sample] == 1:
                samples.add(sample)

        covered_samples.intersection_update(samples)

   
    return covered_samples

def read_combinations(combination_file):
    combinations = []
    with open(combination_file, "r") as reader:
        for line in reader:
            if line == "\n":
                continue
            tokens = line.split()
            combination = tokens[1:5]
            combinations.append(combination)
    return combinations

def compute_stats(TP, FN, TN, FP):
    sensitivity = float(TP) / float(TP + FN)
    specificity = float(TN) / float(TN + FP)
    
    ci_sen = ci(TP, TP + FN, alpha=0.05, method='beta')
    ci_spec = ci(TN, TN + FP, alpha=0.05, method='beta')
    return sensitivity, specificity, ci_sen, ci_spec  


def aggregate_statistics(tumor_matrix, normal_matrix, combinations, gene2id):
    TP = FP = TN = FN = 0
    covered_tumor_samples = set()
    covered_normal_samples = set()
    num_tumor_samples = len(tumor_matrix[0])
    num_normal_samples = len(normal_matrix[0])

    for combination in combinations:
       covered_tumor_samples.update(get_covered_samples(tumor_matrix, gene2id, combination))
       covered_normal_samples.update(get_covered_samples(normal_matrix, gene2id, combination))
     
    TP = len(covered_tumor_samples)
    FN = num_tumor_samples - TP
    FP = len(covered_normal_samples)
    TN = num_normal_samples - FP
    sen, spec, ci_sen, ci_spec = compute_stats(TP, FN, TN, FP)
    return sen, spec, ci_sen, ci_spec, TP, FN, TN, FP

def get_performance(tumor_file, normal_file, combination_file):
    tumor_matrix, gene2id = loadGeneSampleMatrixTumor(tumor_file)
    normal_matrix = loadGeneSampleMatrixNormal(normal_file, gene2id)
    combinations = read_combinations(combination_file)
    return aggregate_statistics(tumor_matrix, normal_matrix, combinations, gene2id)

def gather_performance():
    trTP = trFN = trFP = trTN = 0
    tsTP = tsFN = tsFP = tsTN = 0
    cancer_types = []
    with open("genecount.txt", "r") as reader:
        for line in reader:
            cancer = line.split(",")[0]
            if cancer == "KICH":
                continue
            cancer_types.append(cancer)

    writer.write("Cancer,Training Sensitivity,Train-sen-ci-low,Train-sen-ci-high,Training Specificity,Train-spec-ci-low,Train-spec-ci-high,Test Sensitivity,Test-sen-ci-low,Test-sen-ci-high,Test Specificity,Test-spec-ci-low,Test-spec-ci-high\n")
    for cancer in cancer_types:
        data_dir =  "/gpfs/alpine/world-shared/gen011/sajal/tcga-data/maf2dat-moderate/"
        tumor_test_file = data_dir + "/" + cancer + ".maf2dat.matrix.out.test"
        normal_test_file = data_dir + "/manifest_normal_normal.txt.test.txt.geneSampleList" 


        tumor_train_file = data_dir + "/" + cancer + ".maf2dat.matrix.out.training"
        normal_train_file = data_dir + "/manifest_normal_normal.txt.training.txt.geneSampleList"
        combination_file = "../results/100nodes/" + cancer + "-combinations"

        train_sen, train_spec, train_sen_ci, train_spec_ci, trtp, trfn, trtn, trfp = get_performance(tumor_train_file, normal_train_file, combination_file)
        test_sen, test_spec, test_sen_ci, test_spec_ci, tstp, tsfn, tstn, tsfp = get_performance(tumor_test_file, normal_test_file, combination_file)
        test_sen_low, test_sen_high = test_sen_ci
        test_spec_low, test_spec_high = test_spec_ci
        train_sen_low, train_sen_high = train_sen_ci
        train_spec_low, train_spec_high = train_spec_ci
        
        trTP += trtp
        trFN += trfn
        trTN += trtn
        trFP += trfp
        tsTP += tstp
        tsFN += tsfn
        tsTN += tstn
        tsFP += tsfp
       


        writer.write(",".join([cancer, str(train_sen), str(train_sen_low), str(train_sen_high), str(train_spec), str(train_spec_low), str(train_spec_high), str(test_sen), str(test_sen_low), str(test_sen_high), str(test_spec), str(test_spec_low), str(test_spec_high)]) + "\n")
        print(",".join([cancer, str(train_sen), str(train_spec), str(test_sen), str(test_spec)]) + "\n")
    ci_sen_tr = ci(trTP, trTP + trFN, alpha=0.05, method='beta')
    ci_spec_tr = ci(trTN, trTN + trFP, alpha=0.05, method='beta') 
    ci_sen_ts = ci(tsTP, tsTP + tsFN, alpha=0.05, method='beta')
    ci_spec_ts = ci(tsTN, tsTN + tsFP, alpha=0.05, method='beta')
    sen_tr, spec_tr, ci_sen_tr, ci_spec_tr = compute_stats(trTP, trFN, trTN, trFP)
    sen_ts, spec_ts, ci_sen_ts, ci_spec_ts = compute_stats(tsTP, tsFN, tsTN, tsFP)

    print("Train 95% CI Sensitivity: ", ci_sen_tr, sen_tr)
    print("Train 95% CI Specificity: ", ci_spec_tr, spec_tr)

    print("Test 95% CI Sensitivity: ", ci_sen_ts, sen_ts)
    print("Test 95% CI Specificity: ", ci_spec_ts, spec_ts)
    writer.close()

gather_performance()
