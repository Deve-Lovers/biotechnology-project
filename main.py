import random
import xml.etree.ElementTree as ET


S = "S"
W = "W"
R = "R"
Y = "Y"


probe_size = None
unfinished = []
results_exact = []
results_greedy = []
debug = False
all_sequences = {}
rejected = 0


def weighted_choice(array):
    elements = [item[0] for item in array]
    weights = [item[1] for item in array]
    return random.choices(elements, weights, k=1)[0]


class Probabilities:
    def __init__(self) -> None:
        self.y_a_weight = 0
        self.y_g_weight = 0
        self.y_t_weight = 0
        self.y_c_weight = 0
        self.r_a_weight = 0
        self.r_g_weight = 0
        self.r_t_weight = 0
        self.r_c_weight = 0

    def get_by_chars(self, first_position, last_position) -> int:
        if first_position == "Y" and last_position == "A":
            return self.y_a_weight
        elif first_position == "Y" and last_position == "G":
            return self.y_g_weight
        elif first_position == "Y" and last_position == "T":
            return self.y_t_weight
        elif first_position == "Y" and last_position == "C":
            return self.y_c_weight
        elif first_position == "R" and last_position == "A":
            return self.r_a_weight
        elif first_position == "R" and last_position == "G":
            return self.r_g_weight
        elif first_position == "R" and last_position == "T":
            return self.r_t_weight
        elif first_position == "R" and last_position == "C":
            return self.r_c_weight


def read_xml(file_path):
    try:
        return ET.parse(file_path).getroot()
    except ET.ParseError as e:
        print(f"Error parsing XML: {e}")
        return None


def process_dna_data(dna_data):
    try:
        key = dna_data.attrib['key']
        length = int(dna_data.attrib['length'])
        start_sequence = dna_data.attrib['start']
    except KeyError as e:
        print(f"Missing attribute in DNA data: {e}")
        raise Exception("-")

    probes = dna_data.findall('probe')
    processed_data = []
    for probe in probes:
        pattern = probe.attrib.get('pattern', None)
        cells = [cell.text for cell in probe.findall('cell')]
        processed_data.append((pattern, cells))

    return key, length, start_sequence, processed_data


def create_graph():
    graph_RY, nucleotides = {}, ['A', 'C', 'G', 'T']
    for n in nucleotides:
        graph_RY[n] = []
    return graph_RY


def check_match_ry(pattern, ry_probe):
    conversion_dict_ry = {'A': 'R', 'G': 'R', 'C': 'Y', 'T': 'Y'}
    ry_pattern = ''.join(conversion_dict_ry[n] for n in pattern)
    return ry_pattern == ry_probe[:-1]


def check_match(pattern, ws_probe):
    conversion_dict_ws = {'A': 'W', 'T': 'W', 'C': 'S', 'G': 'S'}
    ws_pattern = ''.join(conversion_dict_ws[n] for n in pattern)
    if ws_pattern == ws_probe[:-1]:
        return True, ws_probe[-1]
    else:
        return False, None
    

def rec_process_graph_greedy(sequence, length, ws_probes, ry_probes, graph_RY, probabilities):
    if len(sequence) == length:
        global results_greedy
        results_greedy.append(sequence)
    else:
        temporary = sequence[-probe_size:]
        optional_nexts = []
        for option in ws_probes:
            next_nuc_check, next_nuc = check_match(temporary, option)
            if next_nuc_check:
                optional_nexts.append(next_nuc)
        mapper = []
        for optional_next in optional_nexts:
            for temp_ry_pattern in graph_RY[optional_next]:
                if check_match_ry(temporary, temp_ry_pattern):
                    weight = probabilities.get_by_chars(temp_ry_pattern[0], optional_next)

                    mapper.append([optional_next, weight])
                    break

        target = weighted_choice(mapper)
        nexts = [target]
        for n in nexts:
            if debug:
                print([nexts, n, sequence, temporary])
            rec_process_graph_greedy(
                f"{sequence}{n}",
                length,
                ws_probes,
                ry_probes=ry_probes,
                graph_RY=graph_RY,
                probabilities=probabilities,
            )

        if nexts == []:
            global unfinished
            unfinished.append(sequence)


def rec_process_graph_exact(sequence, length, ws_probes, ry_probes, graph_RY):
    if len(sequence) == length:
        global results_exact
        results_exact.append(sequence)
    else:
        temporary = sequence[-probe_size:]
        optional_nexts = []
        for option in ws_probes:
            next_nuc_check, next_nuc = check_match(temporary, option)
            if next_nuc_check:
                optional_nexts.append(next_nuc)
        nexts = []
        for optional_next in optional_nexts:
            for temp_ry_pattern in graph_RY[optional_next]:
                if check_match_ry(temporary, temp_ry_pattern):
                    global all_sequences
                    all_sequences[temp_ry_pattern] = True
                    nexts.append(optional_next)
                    break

        global rejected
        if len(optional_nexts) != len(nexts):
            missing = len(optional_nexts) - len(nexts)
            rejected += missing

        for n in nexts:
            if debug:
                print([nexts, n, sequence, temporary])
            rec_process_graph_exact(
                f"{sequence}{n}",
                length,
                ws_probes,
                ry_probes=ry_probes,
                graph_RY=graph_RY,
            )

        if nexts == []:
            global unfinished
            unfinished.append(sequence)


def process_input_data(dna_data, mode: str):
    _, length, start_sequence, probes = process_dna_data(dna_data)

    graph_RY = create_graph()
    ws_probes, ry_probes = probes[0][1], probes[1][1]

    # set size based on single probe
    global probe_size
    probe_size = len(ws_probes[0])-1

    # process graph_RY
    global all_sequences
    for probe in ry_probes:
        graph_RY[probe[-1]].append(probe)
        all_sequences[probe] = False
    for probe in ws_probes:
        all_sequences[probe] = False

    if mode == "exact":
        rec_process_graph_exact(
            sequence=start_sequence,
            length=length,
            ws_probes=ws_probes,
            ry_probes=ry_probes,
            graph_RY=graph_RY,
        )
    else:
        # probabilities
        probabilities = Probabilities()
        for ry_probe in ry_probes:
            if ry_probe[0] == "Y" and ry_probe[-1] == "A":
                probabilities.y_a_weight += 1
            elif ry_probe[0] == "Y" and ry_probe[-1] == "G":
                probabilities.y_g_weight += 1
            elif ry_probe[0] == "Y" and ry_probe[-1] == "T":
                probabilities.y_t_weight += 1
            elif ry_probe[0] == "Y" and ry_probe[-1] == "C":
                probabilities.y_c_weight += 1
            elif ry_probe[0] == "R" and ry_probe[-1] == "A":
                probabilities.r_a_weight += 1
            elif ry_probe[0] == "R" and ry_probe[-1] == "G":
                probabilities.r_g_weight += 1
            elif ry_probe[0] == "R" and ry_probe[-1] == "T":
                probabilities.r_t_weight += 1
            elif ry_probe[0] == "R" and ry_probe[-1] == "C":
                probabilities.r_c_weight += 1
        rec_process_graph_greedy(
            sequence=start_sequence,
            length=length,
            ws_probes=ws_probes,
            ry_probes=ry_probes,
            graph_RY=graph_RY,
            probabilities=probabilities,
        )

    if mode == "exact":
        return results_exact
    return results_greedy    


def main():
    file_path = 'bio2.xml'
    root = read_xml(file_path)
    if root is None:
        print("Failed to read XML file.")
        return
    results = process_input_data(root, mode="exact")
    if debug:
        print("==== unfinished ====")
        if len(unfinished) == 0:
            print("-")
        for i, item in enumerate(unfinished):
            print(i)
            print(item)
    print("==== results ====")
    print(results)
    global rejected
    global all_sequences
    count_true = 0
    count_false = 0
    for _, v in all_sequences.items():
        if v is True:
            count_true += 1
        else:
            count_false += 1
    print(f"Count True: {count_true}")
    print(f"Count False: {count_false}")
    print(f"Rejected {rejected}")
    rejected = 0
    all_sequences = {}
    results = process_input_data(root, mode="greedy")
    if debug:
        print("==== unfinished ====")
        if len(unfinished) == 0:
            print("-")
        for i, item in enumerate(unfinished):
            print(i)
            print(item)
    print("==== results ====")
    print(results)
    count_true = 0
    count_false = 0
    for _, v in all_sequences.items():
        if v is True:
            count_true += 1
        else:
            count_false += 1
    print(f"Count True: {count_true}")
    print(f"Count False: {count_false}")


if __name__ == '__main__':
    main()
