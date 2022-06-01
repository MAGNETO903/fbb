# Что делает эта программа?
# ВАЖНО: РАБОТАЕТ С C6
# Фактически она преобразует файл main_6_2_output.json таким образом,
# чтобы там содержались не B16/B26, а C6 т.е. на выходе мы получим файл такого вида
# у нас есть словарь, где ключи - ncbi ID
# а значения 11 отрезков C6 (отрезков геномов)
# пример:
# {
# <ncbi_id>: [
#               <C6_1>,
#               <C6_2>,
#               <C6_3>,
#               <C6_4>,
#               <C6_5>,
#               <C6_6>,
#               <C6_7>,
#               <C6_8>,
#               <C6_9>,
#               <C6_10>,
#               <C6_11>,
#           ]
# }

import json

# читаем выходной файл main_6_2.py
with open("main_6_2_output.json") as jsonFile:
   all_B26 = json.load(jsonFile)

all_C6 = {}

for i in all_B26:
    cur_C6 = []
    for j in range(len(all_B26[i])):
        if j == 0:
            B16 = all_B26[i][j]
            C6 = B16[:6]
            cur_C6.append(C6)
        else:
            B26 = all_B26[i][j]
            C6 = B26[10:16]
            cur_C6.append(C6)

    all_C6[i] = cur_C6

with open('main_6_6_output.json', 'w') as f:
    f.write(json.dumps(all_C6))
