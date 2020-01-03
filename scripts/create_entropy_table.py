
# импортируем функции из entropy_functions.py в среду, в которой работаем сейчас
import csv
import random
from multiprocessing import Pool
from Bio import SeqIO
from entropy_functions import *

# считывание сиквенса
path_to_seq = "/home/er/Work/entropy-dimers/Y_hg19.fa"  # путь к анализируемому fasta-сиквенсу
path_to_str = "/home/er/Work/entropy-dimers/lobSTR_ystr_hg19.bed" # путь к координатам str
path_to_exon = "/home/er/Work/entropy-dimers/ygenes-exons.txt" # путь к координатам экзонов

# читаем сиквенс с помощью SeqIO
with open(path_to_seq) as handle:
    data = list(SeqIO.parse(handle, "fasta"))
    seq = data[0].seq
    seq = str(seq)

# создадим "базу"
step = 20
threshold = 0.5
sequence = seq
# делим сиквенс на "окна"
segments = {(i, i + step): sequence[i: i + step] for i in range(0, len(sequence), step)}

# используем мультипоточность для ускорения создания базы
pool = Pool()
# считаем энтропии
entropies_raw = pool.map(EntropyRecord.calc_seq_entropy, list(segments.values()))
entropies = dict(zip(list(segments.keys()), entropies_raw))  # энтропии с указанием позиций - "база"

pool.close()
pool.join()

# считываем данные по str в list
str_bed = []
with open(path_to_str) as f:
    for line in f:
        str_bed.append(line.strip().split())

# считываем данные по экзонам в list
exon_file = []
with open(path_to_exon) as f:
    for line in f:
        exon_file.append(line.strip().split())
del exon_file[0] # удаляем строчку с названиями и расшифровкой

str_data = []
# вычисляем энтропии для str
for region_info in str_bed:
    region_name = region_info[5]
    region_type = "STR"
    region_start = region_info[1]
    region_end = region_info[2]
    region_length = int(region_end) - int(region_start)
    if region_length >= 50:
        region_entropy = calculate_entropy(entropies, int(region_start), int(region_end))[1]
        str_data.append([region_name, region_type, region_start, region_end, region_length, region_entropy])

exon_file_sampled = random.sample(exon_file, 10)

exon_data = []
# вычисляем энтропии для экзонов
for region_info in exon_file_sampled:
    region_name = region_info[0]
    region_type = "Protein-exon"
    region_start_positions = list(filter(None, region_info[8].split(',')))
    region_end_positions = list(filter(None, region_info[9].split(',')))
    for i in range(0, len(region_start_positions)):
        region_start = region_start_positions[i]
        region_end = region_end_positions[i]
        region_length = int(region_end) - int(region_start)
        if region_length >= 50:
            region_entropy = calculate_entropy(entropies, int(region_start), int(region_end))[1]
            exon_data.append([region_name, region_type, region_start, region_end, region_length, region_entropy])

exon_data_sampled = random.sample(exon_data, 33)

# TODO? sequence, подаваемая на вход, должна быть уже "обрезанной"
# TODO? или посчитать энтропии заранее, функцией вычислять только median values

# запишем в файл
with open('str_exon_entropy.csv', 'w', newline='') as csvfile:
    table = csv.writer(csvfile, delimiter=' ')
    table.writerow("#name type start end length entropy")
    for entry in str_data:
        table.writerow(entry)
    for entry in exon_data_sampled:
        table.writerow(entry)



# ----------------------------
# использованные функции

def calculate_entropy(entropies, pos_start, pos_end):

    # считаем средние значения в диапазоне pos
    pos = (int(pos_start), int(pos_end))
    median_values = EntropyRecord.calc_pos_entropy(entropies, pos)

    # entropies = EntropyRecord.calc_entropy_over_step(seq, step)
    # pos = (pos_start, pos_end)
    # median_values = EntropyRecord.calc_pos_entropy(entropies, pos)

    return(median_values)
