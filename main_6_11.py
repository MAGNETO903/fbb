# Что же делает эта программа?
# она смотрим на группы B36/B66 и выдаёт строку по каждой из 11 групп в которой содержится число отличий
# на позиции от самого массового B36/B66

import json

# читаем выходной файл main_6_3.py (группы B66)
with open("main_6_3_output.json") as jsonFile:
   B66_groups = json.load(jsonFile)



def key_func(el):
    return el[1]

for i in range(len(B66_groups)):
    print('Группа ', i)
    all_B66 = []
    for B66 in B66_groups[i]:
        all_B66.append([B66, B66_groups[i][B66]['count']])
    all_B66.sort(key=key_func, reverse=True)
    diffs = []

    biggest_group = all_B66[0]
    #print(biggest_group)
    for j in range(len(biggest_group[0])):
        diffs.append(0)
    for j in range(len(all_B66)):
        if j != 0:
            for k in range(len(all_B66[j][0])):
                if biggest_group[0][k] != all_B66[j][0][k]:
                    diffs[k] +=  all_B66[j][1]

    print(diffs)



