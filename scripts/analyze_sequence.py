# скрипт для создания базы энтропий с заданной шириной окна и подсчёта средней энтропии в заданном диапазоне
# использование:
# analyze_sequence.py
# [путь к сиквенсу] [ширина окна] [максимальная доля N]
# [начальная позиция для подсчёта средней энтропии] [конечная позиция]

# импортируем функции из entropy_functions.py в среду, в которой работаем сейчас
import sys
import numpy as np
from multiprocessing import Pool
from Bio import SeqIO
from entropy_functions import *

# проверка на количество аргументов
if len(sys.argv) < 6:
    print(
        "Слишком мало аргументов, нужно 5: "
        "[путь к сиквенсу] [ширина окна] [доля N] [начальная позиция для подсчёта средней энтропии] [конечная позиция]")
    exit()

# считывание сиквенса
path_to_seq = sys.argv[1]  # путь к анализируемому fasta-сиквенсу

# читаем сиквенс с помощью SeqIO
with open(path_to_seq) as handle:
    data = list(SeqIO.parse(handle, "fasta"))
    sequence = data[0].seq

# пример использования функций

# создаём базу
step = int(sys.argv[2])  # ширина "окна"

# укажем минимальную долю не-N для подсчёта энтропии
threshold = float(sys.argv[3])

# переведем сиквенс из типа Seq в string для удобства
sequence = str(sequence)

# делим сиквенс на "окна"
segments = {(i, i + step): sequence[i: i + step] for i in range(0, len(sequence), step)}

# используем мультипоточность для ускорения создания базы
pool = Pool()
# считаем энтропии
entropies_raw = pool.map(EntropyRecord.calc_seq_entropy, list(segments.values()))
entropies = dict(zip(list(segments.keys()), entropies_raw))  # энтропии с указанием позиций - "база"

pool.close()
pool.join()

# считаем средние значения в диапазоне pos
pos = (int(sys.argv[4]), int(sys.argv[5]))
median_values = EntropyRecord.calc_pos_entropy(entropies, pos)

# создаём "базу" по указанным позициям
entropies_base = EntropyRecord.create_pos_entropy_base(entropies, pos)
# записываем базу в файл в текущей директории
np.save("entropies-base.npy", entropies_base)
# файл восстанавливается через read_dictionary = np.load('entropies-base.npy').item()

# выводим полученные данные
print("по данной последовательности:\n", "cреднее количество пройденных пар в окне:", median_values[0], "\n",
      "среднее значение энтропии в диапазоне", pos, ":", "%.5f" % median_values[1], "\n",
      "файл Numpy, содержащий базу в данном диапазоне, сохранён в текущей директории")
