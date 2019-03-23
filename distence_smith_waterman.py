import numpy as np

def smith_waterman_distance(seq1, seq2, match=3, mismatch=-1, insertion=-0.2, deletion=-1, normalize=1):
    tao = 0
    if len(seq2) > len(seq1):
        tao = len(seq1) / len(seq2)
    else:
        tao = len(seq2) / len(seq1)

    if len(seq2) > len(seq1):
        seq1, seq2 = seq2, seq1  # seq1更长
    # 得分矩阵 mat
    mat = np.zeros((len(seq2) + 1, len(seq1) + 1))
    # 列迭代赋值
    for i in range(1, mat.shape[0]):
        # 行迭代赋值
        for j in range(1, mat.shape[1]):
            mat[i, j] = max(
                0,
                # if previous character matches increase the score by match, else decrease it by mismatch
                mat[i - 1, j - 1] + (match if seq1[j - 1] == seq2[i - 1] else mismatch),
                # one character is missing in seq2, so decrease the score by deletion
                mat[i - 1, j] + deletion,
                # one additional character is in seq2, so decrease the scare by insertion
                mat[i, j - 1] + insertion
            )
    # the maximum of mat is now the score, which is returned raw or normalized (with a range of 0-1)
    return np.max(mat) * tao / (len(seq2) * match) if normalize else np.max(mat)

if __name__ == '__main__':
    seq1 = ["a","b","c","d","e","f"]
    seq2 = ["a","b","c","f","a","c"]
    value = smith_waterman_distance(seq1,seq2)
    print(value)


"""
[
    {
        "0":["a","b","c","d","e","f"]
    },
    {
        "1":["d","b","e","f","a","c"]
    },
    {
        "2":["ab","b","c","r","e","h"]
    },
    {
        "3":["a","b","c","f","a","c"]
    }
]
"""