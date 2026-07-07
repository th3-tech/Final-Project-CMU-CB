import argparse

# propensity_file = 'data/propensities.txt'
# sequence_file = 'data/sequence.fasta'

parser = argparse.ArgumentParser(description='Predict the secondary structure of a protein sequence.')
parser.add_argument('sequence_file', type=str, help='Input filename for sequence in data dir')
parser.add_argument('propensities_file', nargs='?', type=str, help='Input filename for propensities in data dir',
                    default='propensities.txt')

args = parser.parse_args()

propensity_file = 'data/' + args.propensities_file
sequence_file = 'data/' + args.sequence_file


def read_sequence(file):
    file = open(file, 'r')
    aa_sequence = ''
    for line in file:
        if line.startswith('>'):
            pass
        else:
            line = line.replace('\n', '').replace('\r', '')
            aa_sequence += line
    return aa_sequence


# creates dictionary and adds AS (key) + desired colum (value) into it
def read_propensities_into_dic(file, column):
    file = open(file, 'r')
    dic_temp = {}
    word_matrix = []
    word = ''
    for line in file:
        if line.startswith('-'):
            pass
        else:
            word_list = []
            for char in line:
                if char != ' ' and char != '\n' and char != '\r':  # \r for linux
                    word = word + char
                elif (char == ' ' or char == '\n' or char == '\r') and word != '':
                    # last value missing if it's last line of file. fixed by adding another '--' line.
                    word_list.append(word)
                    word = ''
                else:
                    pass
            word_matrix.append(word_list)
    for element in word_matrix:
        dic_temp[element[0]] = float(element[column])
    return dic_temp


def map_dic_to_seq(sequence, dictionary):
    value_list = []
    for aa in sequence:
        value_list.append(dictionary[aa])
    return value_list


def find_extendindex(motif, strseq):
    index_list = []
    stpoint = 0
    position = 0
    while position >= 0:
        position = strseq.find(motif, stpoint, len(strseq))
        if position >= 0:
            index_list.append(position)
        stpoint = position + len(motif)
    return index_list


def find_nucleations(value_list, secondary_structure):
    if secondary_structure == 'H':
        window_size = 6
        amount_for_nuc = 4
        nucleation_threshold = 1.03
    elif secondary_structure == 'E':
        window_size = 5
        amount_for_nuc = 3
        nucleation_threshold = 1.00
    sec_struct_list = ['-'] * len(value_list)
    for aa in range(len(value_list) - window_size):
        window = []
        for win_index in range(aa, aa + window_size):
            window.append(value_list[win_index])
        result = list(filter(lambda x: x >= nucleation_threshold, window))
        if len(result) >= amount_for_nuc:
            for i in range(aa, aa + window_size):
                sec_struct_list[i] = secondary_structure
    stringseq = ''.join(sec_struct_list)
    fwd = secondary_structure * 3 + '-'
    rev = '-' + secondary_structure * 3
    # print('Liste von:', secondary_structure)
    # print('fwd', find_extendindex(fwd, stringseq))
    # print('rev', find_extendindex(rev, stringseq))
    checkseq = stringseq
    loopcheck = True
    while stringseq != checkseq or loopcheck is True:
        fwdindex = find_extendindex(fwd, stringseq)
        revindex = find_extendindex(rev, stringseq)
        checkseq = stringseq
        loopcheck = False
        # Extension
        for i in fwdindex:  # fwd
            fwd_threshold = 0
            for j in range(i, i + 4):
                fwd_threshold = fwd_threshold + value_list[j]
            fwd_threshold = fwd_threshold / 4
            if fwd_threshold > 1.00:
                stringseq = stringseq[:i + 1] + secondary_structure * 4 + stringseq[i + 5:]

        for i in revindex:  # rev
            rev_threshold = 0
            for j in range(i, i + 4):
                rev_threshold = rev_threshold + value_list[j]
            rev_threshold = rev_threshold / 4
            if rev_threshold > 1.00:
                stringseq = stringseq[:i + 1] + secondary_structure * 4 + stringseq[i + 5:]
        sec_struct_list = list(stringseq)
    return sec_struct_list


def predict_turns(turn_value_list):
    turn_struct_list = ['-'] * len(turn_value_list)
    threshold = 7.5e-5
    for aa in range(len(turn_value_list) - 3):
        probability = f_i_value_list[aa] * f_i1_value_list[aa + 1] * f_i2_value_list[aa + 2] * f_i3_value_list[aa + 3]
        over1 = (turn_value_list[aa] + turn_value_list[aa + 1] + turn_value_list[aa + 2] + turn_value_list[aa + 3]) / 4
        if probability > threshold and over1 > 1:
            turn_struct_list[aa] = 'T'
    return turn_struct_list


def compile_sec_struct():
    secondary_structure = ['-'] * len(helix_value_list)
    for aa in range(len(secondary_structure)):
        if helix_nuc_list[aa] == 'H' and turn_nuc_list[aa] == '-' and sheet_nuc_list[aa] == '-':  # only H
            secondary_structure[aa] = 'H'
        elif helix_nuc_list[aa] == '-' and turn_nuc_list[aa] == 'T' and sheet_nuc_list[aa] == '-':  # only T
            secondary_structure[aa] = 'T'
        elif helix_nuc_list[aa] == '-' and turn_nuc_list[aa] == '-' and sheet_nuc_list[aa] == 'E':  # only E
            secondary_structure[aa] = 'E'
        elif helix_nuc_list[aa] == 'H' and turn_nuc_list[aa] == '-' and sheet_nuc_list[aa] == 'E':  # H or E
            if helix_value_list[aa] > sheet_value_list[aa]:
                secondary_structure[aa] = 'H'
            else:
                secondary_structure[aa] = 'E'
        elif helix_nuc_list[aa] == '-' and turn_nuc_list[aa] == '-' and sheet_nuc_list[aa] == '-':  # none -> C
            secondary_structure[aa] = 'C'
        elif helix_nuc_list[aa] == 'H' and turn_nuc_list[aa] == 'T' and sheet_nuc_list[aa] == 'E':  # H ore T or E
            alpha = helix_value_list[aa]
            beta = sheet_value_list[aa]
            gamma = turn_value_list[aa]
            max_val = max(alpha, beta, gamma)
            if max_val == alpha:
                secondary_structure[aa] = 'H'
            elif max_val == beta:
                secondary_structure[aa] = 'E'
            elif max_val == gamma:
                secondary_structure[aa] = 'T'
        elif helix_nuc_list[aa] == 'H' and turn_nuc_list[aa] == 'T' and sheet_nuc_list[aa] == '-':  # H or T
            if helix_value_list[aa] > turn_value_list[aa]:
                secondary_structure[aa] = 'H'
            else:
                secondary_structure[aa] = 'T'
        elif helix_nuc_list[aa] == '-' and turn_nuc_list[aa] == 'T' and sheet_nuc_list[aa] == 'E':  # E or T
            if sheet_value_list[aa] > turn_value_list[aa]:
                secondary_structure[aa] = 'E'
            else:
                secondary_structure[aa] = 'T'
    sec_struct_as_string = ''.join(secondary_structure)
    return sec_struct_as_string


sequence = read_sequence(sequence_file)
print('Sequence:'.ljust(10), sequence)

helix_dic = read_propensities_into_dic(propensity_file, 1)
sheet_dic = read_propensities_into_dic(propensity_file, 2)
turn_dic = read_propensities_into_dic(propensity_file, 3)
# print('Helix: ', helix_dic)
# print('Sheet: ', sheet_dic)
# print('Turn: ', turn_dic)
f_i = read_propensities_into_dic(propensity_file, 4)
f_i1 = read_propensities_into_dic(propensity_file, 5)
f_i2 = read_propensities_into_dic(propensity_file, 6)
f_i3 = read_propensities_into_dic(propensity_file, 7)

helix_value_list = (map_dic_to_seq(sequence, helix_dic))
sheet_value_list = (map_dic_to_seq(sequence, sheet_dic))
turn_value_list = (map_dic_to_seq(sequence, turn_dic))
# print('Helix_value_list: ', helix_value_list)
# print('Sheet_value_list: ', sheet_value_list)
# print('Turn_value_list: ', turn_value_list)
f_i_value_list = map_dic_to_seq(sequence, f_i)
f_i1_value_list = map_dic_to_seq(sequence, f_i1)
f_i2_value_list = map_dic_to_seq(sequence, f_i2)
f_i3_value_list = map_dic_to_seq(sequence, f_i3)

turn_nuc_list = predict_turns(turn_value_list)
helix_nuc_list = find_nucleations(helix_value_list, 'H')
sheet_nuc_list = find_nucleations(sheet_value_list, 'E')

print('Helix:'.ljust(10), ''.join(helix_nuc_list))
print('Sheet:'.ljust(10), ''.join(sheet_nuc_list))
print('Turn:'.ljust(10), ''.join(turn_nuc_list))

final_structure = compile_sec_struct()
print('Structure:'.ljust(10), final_structure)

# ----- TESTS BELOW ----- #

# final_test = 'TEHHHHHHHHHHEEETHTEHHHHHEEEEHHTTCEEEETEHEHEHHHHHHTHEHHHHEHC' #>cath|4_2_0|1oaiA00/561-619

'''
final_test2 = []
for i in final_test:
    final_test2.append(i)
print(final_test2)
print(final_structure)
counter = 0

for i in range(len(final_structure)):
    if final_structure[i] == final_test2[i]:
        counter += 1
'''

# print('% equal to online implementation', counter / len(final_structure) * 100)

'''
# http://www.biogem.org/tool/chou-fasman/index.php
# comparison to webtool with implemented cf
'''
