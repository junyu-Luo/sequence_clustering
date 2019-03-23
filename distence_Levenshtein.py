import Levenshtein
import numpy as np

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
seq1 = ["a","b","c","d","e","f"]
seq2 = ["a","b","c","f","a","c"]
print(Levenshtein.seqratio(seq1, seq2))

# s1 = 'abcd'
# s2 = 'abc'
#
# print(Levenshtein.ratio(s1, s2))



# a = np.array([[1, 2, 3, 4],[5, 6, 7, 8]])
# a[0][1] = 5
# print(a)
# dicts = {'NO_0':['1','XB1','XB2'],'NO_1':['XB3','XB4','XB5'],'NO_2':['XB2','XB1','XB0'],'NO_4':['XB0','XB1','XB2','XB3','XB4']}
#
# a = np.zeros((len(dicts), len(dicts)))
# print(a)