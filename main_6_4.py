# Что делает эта программа?
# она группирует С6 т.е. собирает в одинм массив те у которых С в заданной группе одинаковый
# программа вернёт файл следующего вида:
# это массив из 11 элементов, каждый их элементов обозначает номер группы
# элемент имеет следующий вид:
'''

{
    <С6>: {
        "count": <число геномов у которых такой С6
            <id генома 1>,
            <id генома 2>,
            ...
            <id генома N>
        ]
    }
}

'''
# программа работает с выходным файлом main_6_2.py, его структура описана в программе main_6_2.py

import json

# читаем выходной файл main_6_2.py
with open("main_6_2_output.json") as jsonFile:
   all_B26 = json.load(jsonFile)

C6_groups = [
    {},
    {},
    {},
    {},
    {},
    {},
    {},
    {},
    {},
    {},
    {},
]

proccesed_genomes = 0

for i in all_B26:
    proccesed_genomes += 1
    id = i

    for j in range(len(all_B26[i])):
        # мы работаем с B16
        if j == 0:
            B16 = all_B26[i][j]
            C6 = B16[:6]
            #print('C6->', C6)
            if C6 not in C6_groups[j].keys():
                # впервые встретили такой B16
                C6_groups[j][C6] = {
                    "count": 1,
                    "ncbi_ids": [id]
                }
            else:
                #print(1)
                # не впервые
                C6_groups[j][C6]['count'] += 1
                C6_groups[j][C6]['ncbi_ids'].append(id)
        # работаем с B26
        else:
            B26 = all_B26[i][j]
            C6 = B26[10:16]
            #print('C6', C6)
            if C6 not in C6_groups[j].keys():
                # впервые встретили такой B16
                C6_groups[j][C6] = {
                    "count": 1,
                    "ncbi_ids": [id]
                }
            else:

                # не впервые
                C6_groups[j][C6]['count'] += 1
                C6_groups[j][C6]['ncbi_ids'].append(id)

print(proccesed_genomes)

with open('main_6_4_output.json', 'w') as f:
    f.write(json.dumps(C6_groups))






