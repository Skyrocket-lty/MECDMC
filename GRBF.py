import numpy as np


# 读取集成的SM/miRNA的相似性矩阵
SM = np.loadtxt(r'Dateset/dataset1/SM相似性矩阵.txt', dtype=float)
miRNA = np.loadtxt(r'Dateset/dataset1/miRNA相似性矩阵.txt', dtype=float)


#利用高斯径向基函数处理集成的SM/miRNA相似性
def fGRB(A, sigma):
    m, n = A.shape
    B = np.zeros((m, m))

    for i in range(m):
        for j in range(i+1):
            B[i, j] = np.exp(-np.linalg.norm(A[:, i] - A[:, j])**2 / (2 * sigma**2))
            B[j, i] = B[i, j]

    return B



if __name__ == "__main__":
    A_sm = fGRB(SM, 2)
    A_mm = fGRB(miRNA, 2)
    np.savetxt('Dateset/dataset1/SM_GRBF.txt', A_sm, fmt="%0.9f")
    np.savetxt('Dateset/dataset1/miRNA_GRBF.txt', A_mm, fmt="%0.9f")
.