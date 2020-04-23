WINDOW_SIZE = 50
OVERLAP = 0
import sys
import numpy as np
def window_Seq(length):
    windows_ID = []
    i = 0
    while i < length - 1:
        if (i + WINDOW_SIZE - 1) <= length - 1:
            windows_ID.append((i, i + WINDOW_SIZE - 1))
        i = i + WINDOW_SIZE - OVERLAP
    return windows_ID
def sliding_count(genotype):
    count = []
    length = len(genotype)
    i=0
    while i < length:
        info = (genotype[i],genotype[i+1])
        if "0" in info:
            count.append(3)
        else:
            count.append(info.count("1"))
        i = i + 2
    return count
def distance_count(tped_file,indfile_file):
    tped = np.loadtxt(tped_file, dtype=str)
    indfile = np.loadtxt(indfile_file, dtype=str)
    indname = indfile[:, 0]
    length = np.shape(tped)[0]
    ind_number = int((np.shape(tped)[1] - 4) / 2)
    print(ind_number)
    windows_ID = window_Seq(length)
    distance_file = open(tped_file+".dist", "w")
    window_file = open(tped_file+".window", "w")
    for window in windows_ID:
        sub_tped = tped[window[0]:window[1] + 1, :]
        chr = sub_tped[0, 0]
        start = sub_tped[0, 3]
        end = sub_tped[-1, 3]
        geno_count_matrix = []
        geno_tped = sub_tped[:, 4:]
        for genos in geno_tped:
            count = sliding_count(genos)
            geno_count_matrix.append(count)
        geno_count_matrix = np.array(geno_count_matrix)
        distance_matrix = np.zeros(shape=(ind_number, ind_number))
        for m in range(ind_number):
            for n in range(ind_number):
                if m > n:
                    ind1 = geno_count_matrix[:, m]
                    ind2 = geno_count_matrix[:, n]
                    diff = 0
                    miss = 0
                    for i in range(WINDOW_SIZE):
                        geno_diff = (ind1[i], ind2[i])
                        if 3 in geno_diff:
                            miss = miss + 1
                        else:
                            diff = diff + abs(geno_diff[0] - geno_diff[1])
                    if (miss / WINDOW_SIZE) < 0.5:
                        diff = diff + (diff / (WINDOW_SIZE - miss)) * miss
                        distance_matrix[m, n] = distance_matrix[n, m] = diff
                    else:
                        distance_matrix[m, n] = distance_matrix[n, m] = -1
        if -1 in distance_matrix:
            pass
        else:
            window_file.write(chr + "\t")
            window_file.write(start + "\t")
            window_file.write(end + "\t")
            window_file.write("\n")
            for i in range(ind_number):
                distance_file.write(indname[i] + "\t")
                for j in range(ind_number):
                    distance_file.write(str(distance_matrix[i, j]) + "\t")
                    if j == ind_number - 1:
                        distance_file.write("\n")
    window_file.close()
    distance_file.close()
if __name__ == "__main__":
    tped = sys.argv[1]
    fam = sys.argv[2]
    distance_count(tped,fam)
