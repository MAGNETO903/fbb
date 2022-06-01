# Что делает эта программа?
# строит гистограммы по файлу main_6_3_output.json т.е. сгруппированным по B16/B26
# в гистограмму входят ВСЕ группы БЕЗ УЧЁТА ВЕСА ЭТОЙ ГРУППЫ

import json

import matplotlib.pyplot as plt
import numpy as np

with open("main_6_3_output.json") as jsonFile:
   B66_groups = json.load(jsonFile)

data = [
    [],
    [],
    [],
    [],
    [],
    [],
    [],
    [],
    [],
    [],
    []
]

for i in range(len(B66_groups)):
    for j in B66_groups[i]:
        data[i].append(j)

# теперь преобразуем в частоты всё это (без учета веса каждой группы!)
data_freqs = []

for i in range(11):
    data_freqs.append([])
    if i == 0:
        for j in range(36):
            data_freqs[i].append({
                "A": 0,
                "C": 0,
                "T": 0,
                "G": 0
            })
    else:
        for j in range(66):
            data_freqs[i].append({
                "A": 0,
                "C": 0,
                "T": 0,
                "G": 0
            })

for i in range(len(data)):
    for a66 in data[i]:
        for j in range(len(a66)):
            data_freqs[i][j][a66[j]] += 1

# имея частоты построим гистограммы

def get_hist(segment_num, freqs, NAME):
    data_1 = freqs
    x = [[]]
    freq_A = [[]]
    freq_C = [[]]
    freq_T = [[]]
    freq_G = [[]]
    freq_MAX = [[]]

    for i in range(len(data_1)):

        count = 0
        x.append([])
        freq_A.append([])
        freq_C.append([])
        freq_T.append([])
        freq_G.append([])
        #freq_GAP.append([])
        freq_MAX.append([])
        for j in range(len(data_1[i])):
            #if (max(data_1[i][j]['A'], data_1[i][j]['C'], data_1[i][j]['T'], data_1[i][j]['G']) != 0):

                x[-1].append(count+1)
                freq_A[-1].append(data_1[i][j]['A'])
                freq_C[-1].append(data_1[i][j]['C'])
                freq_T[-1].append(data_1[i][j]['T'])
                freq_G[-1].append(data_1[i][j]['G'])
                #freq_MAX[-1].append(max(data_1[i][j]['A'], data_1[i][j]['C'], data_1[i][j]['T'], data_1[i][j]['G']))
                #print(data_1[i][j]['A'] + data_1[i][j]['C'] +data_1[i][j]['T'] + data_1[i][j]['G'])
                count += 1


    plt.figure(figsize=(16,6))

    '''
    if segment_num == 1:
        plt.title('Встречаемость нуклеотида группы ' + freqs['CS'] + ' на позиции в лидер-отрезке')
    else:
        plt.title('Встречаемость нуклеотида группы  ' + freqs['CS'] + ' на позиции в отрезке №' + str(segment_num-1))
    '''

    plt.bar(x[segment_num], freq_A[segment_num], label = 'A')
    plt.bar(x[segment_num], freq_C[segment_num], label = 'C', bottom=np.array(freq_A[segment_num]))
    plt.bar(x[segment_num], freq_T[segment_num], label = 'T', bottom=(np.array(freq_A[segment_num]) + np.array(freq_C[segment_num])) )
    plt.bar(x[segment_num], freq_G[segment_num], label = 'G', bottom=(np.array(freq_A[segment_num]) + np.array(freq_C[segment_num]) + np.array(freq_T[segment_num])))
    #plt.bar(x[segment_num], freq_GAP[segment_num], label = '-', bottom=(np.array(freq_A[segment_num]) + np.array(freq_C[segment_num]) + np.array(freq_T[segment_num]) + np.array(freq_G[segment_num])))

    # zip joins x and y coordinates in pairs
    for x,y in zip(x[segment_num],freq_MAX[segment_num]):

        label = str(round(y*100)/100)

        plt.annotate(label, # this is the text
                     (x,0), # these are the coordinates to position the label
                     textcoords="offset points", # how to position the text
                     xytext=(0,10), # distance from text to points (x,y)
                     ha='center',
                     size = 7
        ) # horizontal alignment can be left, right or center

    plt.legend(loc='upper right', bbox_to_anchor=(1, 1))


    plt.xticks(np.linspace(1,len(freq_A[segment_num]), len(freq_A[segment_num])))
    plt.show()
    #plt.savefig(NAME)

get_hist(1, data_freqs, 'none')