import sys
import matplotlib.pyplot as plt

# Getting (Units vs frequency) data for all repeats
# Dict with each repeat as key and a list of frequency values
# The first value in the list corresponds to the frequency value at 1 repeat unit
repeats_file = sys.argv[1] #Input repeats file (PERF output)

repeats_data = {} # Frequency data at different unit lengths for a Repeat class
with open(repeats_file) as fh:
    for line in fh:
        line = line.strip().split("\t")
        rep_class = line[3]
        units = int(line[6])
        try:
            try:
                repeats_data[rep_class][units - 1] += 1
            except IndexError:
                repeats_data[rep_class] = repeats_data[rep_class] + [0]*(units - len(repeats_data[rep_class]))
                repeats_data[rep_class][units - 1] += 1
        except KeyError:
            repeats_data[rep_class] = [0]*units
            repeats_data[rep_class][units - 1] += 1

def bell_curve(repeats_data): # Function to detect bell curve
    output = {}
    for rep in sorted(repeats_data.keys()):
        freq_data = repeats_data[rep]
        min_unit = 12//len(rep) # Starting repeat unit to min repeat length
        pf = freq_data[min_unit - 1] # Initialise with freq at min unit
        check = 0 # Check for increase in frequency
        span_freqs = [] # Storing freqs in the span
        for u, f in enumerate(freq_data[min_unit:]):
            u = u + min_unit + 1 # Actual repeat unit number
            if f >= pf:
                if check == 0 and pf > 10: # register min unit, min freq if first increase
                    span_freqs = [pf, f] #storing freqs from the prev unit
                    peak_start, minf, check = [u - 1, pf, 1] # Peak start unit is u-1
                elif check == 1: # storing freq if peak already detected
                    span_freqs.append(f)
            else:
                if check == 1:
                    span_freqs.append(f) # store the least
                    if f < minf and peak_start > min_unit: # if in span and freq is found below min freq
                        # store to output
                        if rep not in output: 
                            output[rep] = [[peak_start, u, u - peak_start + 1, span_freqs]]
                        else:
                            output[rep].append([peak_start, u, u - peak_start + 1, span_freqs])
                        check = 0
            pf = f #updating previous frequency
    return output

output = bell_curve(repeats_data)
for rep in output:
    for i, o in enumerate(output[rep]):
        if i > 0:
            if o[0] == pstop:
                output[rep][i] = [pstart, o[1], o[1] - pstart + 1, pfreqs[:-1] + o[3]]
                output[rep][i-1] = 0
        pstop = o[1]
        pstart = o[0]
        pfreqs = o[3]

print('repeat', 'start', 'stop', 'span', 'amplitude', 'norm_amplitude', 'sum_of_freqs_span', 'percent_of_freqs_in_span', 'normalised_percent', sep='\t')
for rep in output:
    for o in output[rep]:
        if isinstance(o, list) and o[2] >= 4: # if span of the peak greater that three repeat units
            freqs = o[3]
            diff = max(o[3]) - o[3][0]
            norm_diff = round((diff/sum(repeats_data[rep]))*1000, 2)
            span_sum = sum(freqs[1:-1])
            percent_span_sum = round((span_sum/sum(repeats_data[rep]))*100, 2)
            norm_percent_span_sum = round(percent_span_sum/(o[2]-2), 2)
            if diff > 0:
                print(rep, o[0], o[1], o[2], diff, norm_diff, span_sum, percent_span_sum, norm_percent_span_sum, sep="\t")