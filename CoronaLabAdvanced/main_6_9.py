# Что делает эта программа?
# она работает с файлом main_6_3_output.json
# т.е. с группами B36/B66
# выводит в столбик все B66 в каждой из 11 групп и рядом подписывает количество геномов, у которых такой B66

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
    for j in range(len(all_B66)):
        print('> group'+str(i) + ' ' +str(round((all_B66[j][1]/359669)*100, 5))+'% ('+str(all_B66[j][1])+')')
        print(all_B66[j][0])


