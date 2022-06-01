from Bio import SeqIO
from Bio.Seq import Seq

import pandas as pd

import re

from IPython.display import display, HTML


def extract_info(file_name):
    # геном (строка)
    genome_str = ""
    # индентификатор генома (штамм коронавируса)
    genome_id = ""

    # название вируса
    virus_name = ""

    file1 = file_name

    # считыванием генома и названия его штамма
    for record in SeqIO.parse(file1, "genbank"):
        genome_str = str(record.seq)
        genome_id = record.id
        virus_name = record.description

    codes = []

    recs = [rec for rec in SeqIO.parse(file1, "genbank")]
    for rec in recs:
        feats = [feat for feat in rec.features if feat.type == "CDS"]
        for feat in feats:
            # print(feat)
            found = False
            try:
                name = feat.qualifiers['gene'][0]
                found = True
            except KeyError:
                found = False

            if found == False:
                try:
                    name = feat.qualifiers['note'][0]
                    found = True
                except KeyError:
                    found = False
                    # print(feat)

            if found == False:
                try:
                    name = feat.qualifiers['product'][0]
                    found = True
                except KeyError:
                    found = False
                    print(feat)

            codes.append([
                int(feat.location.start),
                int(feat.location.end),
                name
            ])

    return {
        "markup": codes,
        "genome": genome_str,
        "id": genome_id,
        "name": virus_name
    }


def word_sort(el):
    return el['matches_num']


def match_sort(el):
    return el['global_pos']


# лимит длины апстрима
U_LIMIT = 10

main_markup = []
# апстримы
u_markup = []
# f-интервалы
f_markup = []

path_to_file = str(input("Введите имя (или путь до файла) вируса (например, SARS): "))  # "HCoV_gb.gb"
# path_to_file = 'SARS'

path_to_file += '.gb'

info = extract_info(path_to_file)

main_markup = info["markup"]
genom = info["genome"]
genome_id = info["id"]
genome_name = info["name"]

# print(info)

leader_markup = [0, int(main_markup[0][0])]

atg_markup = [m.end() for m in re.finditer('ATG', str(genom))]


# print(atg_markup)
def get_nearest_atg(pos):
    best_atg_dist = 100000
    best_atg_pos = -1
    for j in range(len(atg_markup)):
        dist = pos - atg_markup[j]

        if dist > 0 and dist < best_atg_dist:
            best_atg_pos = atg_markup[j]
            best_atg_dist = dist
    return best_atg_pos


def get_consensus(arr_of_genes, cons_len=18):
    consensus = ''
    # print(cons_len)
    for i in range(cons_len):
        # print(i)
        freq = {
            "C": 0,
            "T": 0,
            "A": 0,
            "G": 0
        }
        for j in range(len(arr_of_genes)):
            freq[arr_of_genes[j][i]] += 1

        # print(freq)

        max_val = 0

        for j in ['C', 'T', 'A', 'G']:
            if freq[j] > max_val:
                max_val = freq[j]

        max_val_num = 0

        for j in ['C', 'T', 'A', 'G']:
            if freq[j] == max_val:
                max_val_num += 1
                if max_val_num == 1:
                    consensus += str(j)
                else:
                    consensus += '/' + str(j)

    return consensus


def get_atg_inside(s, e):
    out = []
    for j in range(len(atg_markup)):
        if atg_markup[j] > s and atg_markup[j] + 2 < e:
            out.append(atg_markup[j])

    return out


# найдём все апстримы
for i in range(1, len(main_markup)):
    criterion = ''

    # начало текущего гена
    start = int(main_markup[i][0])

    name = main_markup[i][-1]

    # нашли, или нет
    found = False

    # конец пред. гена
    prev_gen_end = int(main_markup[i - 1][1])
    # начало пред. гена
    prev_gen_start = int(main_markup[i - 1][0])

    # находим ближайший ATG
    nearest_atg = get_nearest_atg(start)

    # если этот ATG между генами
    if prev_gen_end < nearest_atg < start:
        # берем его если проходим по длине
        if start - nearest_atg > U_LIMIT:
            criterion = 'external_atg'
            u_markup.append([nearest_atg, start, criterion])

            found = True

    # если расстояние между генами достаточное
    if found == False:
        if start - prev_gen_end > U_LIMIT:
            criterion = 'prev_gen_stop'
            u_markup.append([prev_gen_end, start, criterion])

            found = True

    # найдем все атг внутри предыдущего гена
    if found == False:
        atg_inside_prev_gen = get_atg_inside(prev_gen_start, prev_gen_end)

        atg_inside_prev_gen.reverse()
        for atg in atg_inside_prev_gen:
            if start - atg > U_LIMIT:
                criterion = 'interior_atg'
                u_markup.append([atg, start, criterion])

                found = True
                break

    if found == False:
        # на крайний случай от начала пред. гена до текущего гена
        criterion = 'prev_gen_start'
        u_markup.append([prev_gen_start, start, criterion])

    # если сильно длинный апстрим
    if u_markup[-1][1] - u_markup[-1][0] > 60:
        criterion = '60'
        # обрезаем его до 60 символов
        u_markup[-1][0] = u_markup[-1][1] - 60
        u_markup[-1][2] = criterion

for i in range(len(u_markup) - 1):
    f_markup.append([u_markup[i][1] + 1, u_markup[i + 1][0] - 1])

f_markup.append([u_markup[-1][1] + 1, len(str(genom)) - 1])

# print("U- и F-stream'ы:")
for i in range(len(u_markup)):
    # print("U {} перед {}".format(u_markup[i], [main_markup[i+1][0], main_markup[i+1][1], main_markup[i+1][-1]]))
    u_markup[i][0] -= 20
    if i > 0:
        f_markup[i - 1][1] -= 20
    if i < len(u_markup):
        pass
        # print("F",f_markup[i])

# алфавит
alphabet = ['A', 'C', 'T', 'G']

# частота слов
freq = {
}
# сами слова
all_possible_words = []
# перебор всех 6-буквенных слов
for s1 in alphabet:
    for s2 in alphabet:
        for s3 in alphabet:
            for s4 in alphabet:
                for s5 in alphabet:
                    for s6 in alphabet:
                        word = s1 + s2 + s3 + s4 + s5 + s6
                        num = len([m.start() for m in re.finditer(word, str(genom))])
                        all_possible_words.append(word)
                        freq[word] = num

words_in_leader = []
words_candidates = []

# у нас есть все слова, теперь отберем те, что есть в leader
for word in all_possible_words:
    if len([m.start() for m in re.finditer(word, str(genom)[leader_markup[0]:leader_markup[1]])]) > 0:
        words_in_leader.append(word)

# и те, что в upstram'ах
goal_num = int(len(u_markup) * 0.4)

special_words = []
for word in words_in_leader:
    matches_num = 0
    strange = False
    u_matches = []
    l_matches = []
    f_matches = []
    for j in range(len(u_markup)):
        matches = [{"num": j + 1, "local_pos": m.start(), "global_pos": m.start() + u_markup[j][0]} for m in
                   re.finditer(word, str(genom)[u_markup[j][0]:u_markup[j][1]])]
        if len(matches) > 1:
            strange = True
            matches_num += 1
        elif len(matches) == 1:
            matches_num += 1
        u_matches += matches

    for j in range(len(f_markup)):
        matches = [{"num": j + 1, "local_pos": m.start(), "global_pos": m.start() + f_markup[j][0]} for m in
                   re.finditer(word, str(genom)[f_markup[j][0]:f_markup[j][1]])]
        f_matches += matches

    l_matches = [{"local_pos": m.start(), "global_pos": m.start() + leader_markup[0]} for m in
                 re.finditer(word, str(genom)[leader_markup[0]:leader_markup[1]])]

    if matches_num >= goal_num:
        if (strange):
            special_words.append(word)
        else:
            words_candidates.append({
                "word": word,
                "matches_num": matches_num,
                "l_matches": l_matches,
                "u_matches": u_matches,
                "f_matches": f_matches
            })

# print("Главные кандитаты:")


words_candidates.sort(key=word_sort, reverse=True)

'''
for i in range(len(words_candidates)):

    print(words_candidates[i]["word"], " - в {}/{}".format(words_candidates[i]["matches_num"], len(u_markup)), "upstream'ах")
    print("встречи в лидере: {}".format(words_candidates[i]["l_matches"]))
    print("встречи в апстримах: {}".format(words_candidates[i]["u_matches"]))
    print("встречи в F-стримах: {}".format(words_candidates[i]["f_matches"]))

    global_matches = [m.start() for m in re.finditer(words_candidates[i]["word"], str(genom))]
    print("Позиции {} во всем геноме (всего {} вхождений):".format(words_candidates[i]['word'], len(global_matches)))
    print(global_matches)
    print("Выравнивания кандитатов")
    l_match_pos = words_candidates[i]['l_matches'][0]['global_pos']
    matches_18_n = []
    matches_18_id = []
    print(words_candidates[i]['l_matches'][0]['global_pos'])
    if l_match_pos > 6:
        print("L  :", genom[l_match_pos-6: l_match_pos+12])
    else:
        print("L  :", genom[l_match_pos - (6-l_match_pos): l_match_pos + 12])

    matches_18_n.append(genom[l_match_pos-6: l_match_pos+12])
    matches_18_id.append('L')

    for j in range(len(words_candidates[i]['u_matches'])):
        match_pos = words_candidates[i]['u_matches'][j]['global_pos']
        upstream_num = words_candidates[i]['u_matches'][j]['num']
        if upstream_num < 10:
            print("U{} : {}".format(upstream_num, genom[match_pos-6:match_pos+12]))
        else:
            print("U{}: {}".format(upstream_num, genom[match_pos - 6:match_pos + 12]))

        matches_18_n.append(genom[match_pos - 6:match_pos + 12])
        matches_18_id.append('U{}'.format(upstream_num))

    for j in range(len(words_candidates[i]['f_matches'])):
        match_pos = words_candidates[i]['f_matches'][j]['global_pos']
        fstream_num = words_candidates[i]['f_matches'][j]['num']
        if fstream_num < 10:
            print("F{} : {}".format(fstream_num, genom[match_pos-6:match_pos+12]))
        else:
            print("F{}: {}".format(fstream_num, genom[match_pos - 6:match_pos + 12]))

        matches_18_n.append(genom[match_pos - 6:match_pos + 12])
        matches_18_id.append('F{}'.format(fstream_num))

    print("консенсус: \n {}".format(get_consensus(matches_18_n)))
    ofile = open(("{}_{}.fasta".format(genome_id, words_candidates[i]["word"])), "w")

    for j in range(len(matches_18_n)):

        ofile.write(">" + matches_18_id[j] + "\n" + matches_18_n[j] + "\n")

    # do not forget to close it
    ofile.close()

print("Специфичные кандидаты, те, которые в одном из upstream'ах встречаются больше одного раза")
print(special_words)

'''
genome_name = genome_name.split(', complete genome')[0]

genomes = {'название': [], 'начало': [], 'конец': []}
markups = {'название': ['L'], 'начало': [leader_markup[0]], 'конец': [leader_markup[1]], 'критерий': ['-']}

for i in range(len(main_markup)):
    genomes['название'].append(main_markup[i][2])
    genomes['начало'].append(main_markup[i][0])
    genomes['конец'].append(main_markup[i][1])

for i in range(len(u_markup)):
    markups['название'].append("U" + str(i + 1))
    markups['начало'].append(u_markup[i][0])
    markups['конец'].append(u_markup[i][1])
    markups['критерий'].append(u_markup[i][2])

    markups['название'].append("F" + str(i + 1))
    markups['начало'].append(f_markup[i][0])
    markups['конец'].append(f_markup[i][1])
    markups['критерий'].append('-')

# print(genome_id, genome_name)
# print("Гены")
# print(pd.DataFrame(data=genomes))
# display(HTML(pd.DataFrame(data=genomes).to_html()))
# print("Разметка")
# display(HTML(pd.DataFrame(data=markups).to_html()))
# print(pd.DataFrame(data=markups))
# print("Находки")
for i in range(len(words_candidates)):
    # print("{} L:{} U:{} F:{}".format(words_candidates[i]["word"], len(words_candidates[i]["l_matches"]), len(words_candidates[i]["u_matches"]), len(words_candidates[i]["f_matches"])))
    # print("список слов:")
    all_matches = []
    all_matches.append({
        "name": "L",

        "global_pos": words_candidates[i]['l_matches'][0]['global_pos']
    })
    for j in range(len(words_candidates[i]['u_matches'])):
        num = words_candidates[i]['u_matches'][j]['num']
        local_pos = words_candidates[i]['u_matches'][j]['local_pos']
        global_pos = words_candidates[i]['u_matches'][j]['global_pos']

        all_matches.append({
            "name": 'U' + str(num),

            "global_pos": global_pos
        })
    # all_matches += words_candidates[i]['u_matches']
    for j in range(len(words_candidates[i]['f_matches'])):
        num = words_candidates[i]['f_matches'][j]['num']
        local_pos = words_candidates[i]['f_matches'][j]['local_pos']
        global_pos = words_candidates[i]['f_matches'][j]['global_pos']

        all_matches.append({
            "name": 'F' + str(num),

            "global_pos": global_pos
        })
    all_matches.sort(key=match_sort)
    matches = {"участок": [], "позиция в геноме": [], "дистанция до позднего гена": []}
    for match in all_matches:
        # найдем дистанция
        dist_post_gene = 1000000000
        for post_gene in main_markup:
            s = post_gene[0]
            dist = s - match['global_pos'] + 6
            if dist > 0 and dist < dist_post_gene:
                dist_post_gene = dist
        if dist_post_gene == 1000000000:
            dist_post_gene = 0
        matches['участок'].append(match['name'])
        matches['дистанция до позднего гена'].append(dist_post_gene)
        matches['позиция в геноме'].append(match['global_pos'])
    # display(HTML(pd.DataFrame(data=matches).to_html()))

    matches_18_n = []
    matches_18_id = []
    l_match_pos = words_candidates[i]['l_matches'][0]['global_pos']

    matches_18_n.append(genom[l_match_pos - 6: l_match_pos + 12])
    matches_18_id.append('L')

    for j in range(len(words_candidates[i]['u_matches'])):
        match_pos = words_candidates[i]['u_matches'][j]['global_pos']
        upstream_num = words_candidates[i]['u_matches'][j]['num']

        matches_18_n.append(genom[match_pos - 6:match_pos + 12])
        matches_18_id.append('U{}'.format(upstream_num))

    for j in range(len(words_candidates[i]['f_matches'])):
        match_pos = words_candidates[i]['f_matches'][j]['global_pos']
        fstream_num = words_candidates[i]['f_matches'][j]['num']

        matches_18_n.append(genom[match_pos - 6:match_pos + 12])
        matches_18_id.append('F{}'.format(fstream_num))

    ofile = open(("{}_{}.fasta".format(genome_id, words_candidates[i]["word"])), "w")

print("Предсказанная CS - {}:".format(words_candidates[0]["word"]))
# print(words_candidates)

# ТАБЛИЦА ДЛЯ ПОЛЬЗОВАТЕЛЕЙ
for_reading = {
    "имя гена": [],
    "координата старта трансляции": [],
    "обнаружен CS": [],
    "расстояние от CS до старта трансляции": [],
    "примечания": []
}

# print(words_candidates)

best_CS_u_matches = words_candidates[0]['u_matches']
best_CS_u_matches.insert(0, words_candidates[0]['l_matches'][0])
best_CS_u_matches[0]['num'] = 0

for i in range(len(main_markup)):
    name = main_markup[i][2]
    start_pos = main_markup[i][0]
    for_reading['имя гена'].append(name)
    for_reading['координата старта трансляции'].append(start_pos)
    exist = False
    CS_match_pos = ' '
    for j in range(len(best_CS_u_matches)):
        num = best_CS_u_matches[j]['num']
        if num == i:
            exist = True
            CS_match_pos = best_CS_u_matches[j]['global_pos']
            break

    if exist:
        for_reading['обнаружен CS'].append('да')
        dist = (start_pos - CS_match_pos) - 6
        for_reading['расстояние от CS до старта трансляции'].append(dist)
    else:
        for_reading['обнаружен CS'].append('нет')
        for_reading['расстояние от CS до старта трансляции'].append(" ")

    for_reading['примечания'].append(' ')

for_reading['имя гена'][0] = 'leader'

for_reading['имя гена'].pop(1)
for_reading['координата старта трансляции'].pop(1)
for_reading['обнаружен CS'].pop(1)
for_reading['расстояние от CS до старта трансляции'].pop(1)
for_reading['примечания'].pop(1)

display(pd.DataFrame(data=for_reading))
# print(main_markup)