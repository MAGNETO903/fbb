# Что делает эта программа?
# она группирует B36/B66 т.е. собирает в одинм массив те у которых B36/B66 в заданной группе одинаковый
# программа вернёт файл следующего вида:
# это массив из 11 элементов, каждый их элементов обозначает номер группы
# элемент имеет следующий вид:
'''

{
    <B36/B66>: {
        "count": <число геномов у которых такой B16/B26
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

B26_groups = [
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

            if B16 not in B26_groups[j].keys():
                # впервые встретили такой B16
                B26_groups[j][B16] = {
                    "count": 1,
                    "ncbi_ids": [id]
                }
            else:
                #print(1)
                # не впервые
                B26_groups[j][B16]['count'] += 1
                B26_groups[j][B16]['ncbi_ids'].append(id)
        # работаем с B26
        else:
            B26 = all_B26[i][j]

            if B26 not in B26_groups[j].keys():
                # впервые встретили такой B16
                B26_groups[j][B26] = {
                    "count": 1,
                    "ncbi_ids": [id]
                }
            else:

                # не впервые
                B26_groups[j][B26]['count'] += 1
                B26_groups[j][B26]['ncbi_ids'].append(id)

print(proccesed_genomes)

with open('main_6_3_output.json', 'w') as f:
    f.write(json.dumps(B26_groups))






