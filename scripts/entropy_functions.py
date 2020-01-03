# функции для посчёта энтропии димеров (на основе метода из пакета lobSTR)
# примеры использования находятся в конце

from math import log2, floor
import numpy as np

# расширенный словарь из lobstr_index.py (учитывает неоднозначности вроде K, M, R, S, W)
nucToNumber = {"A": 0, "a": 0, "C": 1, "c": 1, "G": 2, "g": 2, "T": 3, "t": 3,
               "W": 3, "w": 3, "B": 3, "b": 3, "H": 3, "h": 3,
               "K": 2, "k": 2, "S": 2, "s": 2, "D": 2, "d": 2,
               "M": 0, "m": 0, "R": 0, "r": 0,
               "Y": 1, "y": 1, "V": 1, "v": 1}


# функция, для подсчёта "нормированной" величины -p*log(p)
def minus_p_log_p(value):
    return -1 * value * log2(value)


# функции на основе С++ кода из EntropyDetection.cpp, lobstr_index.py
# названия аналогичных функций и переменных сохранены
class EntropyDetection:
    # класс, содержит функции для подсчёта энтропии

    def entropy_one_window_dinuc(window_nucs):
        # функция, принимает окно, считает его энтропию
        # возвращает конструкцию (число рассмотренных пар, энтропия)

        window_length = len(window_nucs)
        kmer_counts = np.repeat(0, 34)
        pairs = 0

        for i in range(0, window_length - 1):
            nuc1 = window_nucs[i]
            nuc2 = window_nucs[i + 1]
            if nuc1 != 'N' and nuc2 != 'N':
                kmer_counts[nucToNumber[nuc1] + nucToNumber[nuc2] * 10] += 1
                pairs += 1

        entropy = 0
        for i in range(0, len(kmer_counts)):
            count = kmer_counts[i]
            if count != 0:
                p = count / pairs
                entropy += minus_p_log_p(p)

        return pairs, (4 - entropy) / 4


class EntropyRecord:
    # класс, функции для подсчёта энтропии для выявления границы для STR-регионов
    # THRESHOLD - доля баз не-N, которые должны быть в сиквенсе, чтобы мы могли мерить энтропию
    # step - шаг, с которым двигаем окно, в котором считаем энтропию

    def calc_seq_entropy(seq):
        # функция, считает энтропию последовательности

        global threshold
        threshold = 0.5

        if len(seq) < 2:
            return 0, 0

        # проверка условия на на THRESHOLD
        if threshold <= seq.count('N') / len(seq):
            return 0, 0

        return EntropyDetection.entropy_one_window_dinuc(seq)

    def calc_entropy_over_step(seq, step):
        # функция, считает энтропию последовательности по окнам, двигая окно через данный шаг.
        # Записывает энтропии в словарь (базу данных)
        # возвращает словарь с координатами окон и энтропиями
        # (последовательность индексируется с 0)

        entropies = {}

        for i in range(0, len(seq), step):
            start = i
            end = i + step - 1
            if i + step > len(seq):  # если дошли до последнего члена
                end = len(seq)

            window = seq[start:end]
            window_entropy = EntropyRecord.calc_seq_entropy(window)
            entropies[(start, end)] = window_entropy

        return entropies

    def calc_pos_entropy(entropies, pos):
        # функция, считает усреднённые энтропию и число пар по координатам, заданным в pos
        # pos = (start, end)
        # ширина окна = step
        # entropies - база посчитанных энтропий для последовательности seq

        start = pos[0]
        end = pos[1]

        windows_num = 0
        entropy_sum = 0
        pairs_sum = 0

        for coordinates in entropies.keys():
            if (coordinates[0] <= start < coordinates[1]) or (coordinates[0] < end <= coordinates[1])\
                    or ((start < coordinates[0]) and (end > coordinates[1])):
                windows_num += 1
                pairs_sum += entropies[coordinates][0]  # число пар
                entropy_sum += entropies[coordinates][1]  # энтропия

        return floor(pairs_sum / windows_num), entropy_sum / windows_num

    def create_pos_entropy_base(entropies, pos):
        # функция, создаёт базу энтропий по позициям из словаря entropies
        # по сути, slice
        # pos = (start, end)
        # entropies - база посчитанных энтропий для последовательности seq

        start = pos[0]
        end = pos[1]

        entropies_sliced = {}

        for coordinates in entropies.keys():
            if (coordinates[0] <= start < coordinates[1]) or (coordinates[0] < end <= coordinates[1])\
                    or ((start < coordinates[0]) and (end > coordinates[1])):

                entropies_sliced[coordinates] = (entropies[coordinates][0], entropies[coordinates][1])

        return entropies_sliced

# usage
# threshold меняется во вспомогательной функции, т.к. здесь задан через global (начальное значение = 0.5)
# seq - последовательность
# step - ширина окна
# entropies = EntropyRecord.calc_entropy_over_step(seq, step) создаёт базу entropies
# pos = (start, end) - диапазон, по которому усреднять энтропии
# EntropyRecord.calc_pos_entropy(entropies, pos) считает среддние значения в
# диапазоне pos

# test
# 1) EntropyDetection.entropy_one_window_dinuc("ATGC")
# 2) EntropyRecord.calc_seq_entropy("ATtGC", 0.2)
# 3) EntropyRecord.calc_entropy_over_step("ATGCTGATGCT", 5)
# 4) EntropyRecord.calc_pos_entropy(entropies, pos)

