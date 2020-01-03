# импортируем функции из entropy_functions.py в среду, в которой работаем сейчас
from multiprocessing import Pool
from Bio import SeqIO
from entropy_functions import *
from itertools import chain
import csv
import json

# считывание сиквенса
path_to_seq = "/home/er/Work/entropy-dimers/Y_hg19.fa"  # путь к анализируемому fasta-сиквенсу
path_to_str = "/home/er/Work/entropy-dimers/lobSTR_ystr_hg19.bed" # путь к координатам str
path_to_exon = "/home/er/Work/entropy-dimers/ygenes-exons.txt" # путь к координатам экзонов

# читаем сиквенс с помощью SeqIO
with open(path_to_seq) as handle:
    data = list(SeqIO.parse(handle, "fasta"))
    seq = data[0].seq
    seq = str(seq)

threshold = 0.5

# создадим "базу"
step = 15
cut_left = 10000
cut_right = 30000000
sequence = seq[cut_left:cut_right]

entropies_list = []

# посчитаем несколько баз энтропий, каждая последующая со сдвигом по позиции на 1 относительно нуля
for start in range(0, step):
    entropy = create_entropy_base(sequence, step, start, cut_left)
    sorted_entropy = sorted(entropy.items(), key=lambda kv: kv[0]) # сортируем словарь по позициям, превращая его в list
    #entropies_list.append(sorted_entropy)

    # упростим каждый из получившихся листов для последующей более удобной индексации
    newlist = list(chain.from_iterable(sorted_entropy))
    newlist = list(chain.from_iterable(newlist))
    # оставим только значения энтропий окон, т.к. мы примерно представляем, по каким позициям они идут
    newlist = newlist[3::4]
    entropies_list.append(newlist)

position = 15
summary_entropy = 0

entropies_list_backup = entropies_list
# entropies_list = entropies_list_backup

# пока пустой массив, в который будем записывать энтропии по каждой позиции
entropies_positions = np.repeat(0, cut_right-cut_left-step*2)

pos_dict_final = dict()

for i in range(0, step-1):

    entropy_sum = [sum(x)/step for x in zip(*entropies_list)]

    pos_keys = range(cut_left + i, cut_right + i, step)
    pos_dict = dict(zip(pos_keys, entropy_sum))
    pos_dict_final.update(pos_dict)
    # сдвинем один из листов по индексу (~как сдвиг рамки)
    entropies_list[step-1-i] = [entropies_list[step-1-i][-1]] + entropies_list[step-1-i][:-1]

# отсортируем
pos_dict_final = sorted(pos_dict_final.items(), key=lambda kv: kv[0])

# запишем

with open('/home/er/Work/entropy-dimers/positions-entropy.csv', 'w', newline='') as out:
    csv_out=csv.writer(out, delimiter = ' ')
    csv_out.writerow(['position', 'entropy'])
    for row in pos_dict_final:
        csv_out.writerow(row)

def create_entropy_base(sequence, step, start, offset_coordinate):

    # делим сиквенс на "окна"
    segments = {(i + start + offset_coordinate, i + step + start + offset_coordinate) : sequence[i + start:i + step + start] for i in range(0, len(sequence), step)}

    # используем мультипоточность для ускорения создания базы
    pool = Pool()
    # считаем энтропии
    entropies_raw = pool.map(EntropyRecord.calc_seq_entropy, list(segments.values()))
    entropies = dict(zip(list(segments.keys()), entropies_raw))  # энтропии с указанием позиций - "база"

    pool.close()
    pool.join()

    return (entropies)

#sample_entropy = create_entropy_base(sequence, 15, 3, 10000)