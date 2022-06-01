# Что делает эта программа?
# она работает с файлом main_6_6_output.py
# она ищет те геномы, у которых каждый C6 отличается хотя бы на один символ от референса и выводит ID этого генома
# и соотвествующие C6

import json

# референсные C6
ref_C6 = [
    'ACGAAC',
    'ACGAAC',
    'ACGAAC',
    'ACGAAC',
    'ACGAAC',
    'ACGAAC',
    'ACGAAC',
    'AAGAAA',
    'ACGAAC',
    'ACGAAC',
    'ACAATC'
]

# кол-во несовпадающих символов
def diff(seq1, seq2):
    miss = 0
    for i in range(len(seq1)):
        if seq1[i] != seq2[i]:
            miss += 1
    return miss

# читаем выходной файл main_6_2.py
with open("main_6_6_output.json") as jsonFile:
   all_C6 = json.load(jsonFile)



for i in all_C6:
    diffs = 0
    for j in range(len(all_C6[i])):
        cur_ref_C6 = ref_C6[j]
        cur_C6 = all_C6[i][j]

        #print()

        if diff(cur_C6, cur_ref_C6) > 0:
            diffs += 1

    if diffs > 5:
        print(i)
        print(all_C6[i])