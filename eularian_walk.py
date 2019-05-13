import random
from datetime import datetime
from tqdm import tqdm
from termcolor import colored
import timeit
import matplotlib.pyplot as plt
import numpy as np

def perform_walk(d_in, d_out, in_count, out_count) -> str:
    ## https://www.geeksforgeeks.org/hierholzers-algorithm-directed-graph/
    ## https://math.stackexchange.com/questions/1871065/euler-path-for-directed-graph
    a_n = ""
    b_n = ""
    ## find a_n and b_n
    # a) There should be a single vertex in graph which has (indegree+1==outdegree): An
    # b) There should be a single vertex in graph which has (indegree==outdegree+1): Bn

    for node, _ in in_count.items():
        if in_count[node] + 1 == out_count[node]:
            # print(colored("a_n:", 'magenta', attrs=['bold']), node)
            a_n = node
        elif in_count[node] == out_count[node] + 1:
            # print(colored("b_n:", 'magenta', attrs=['bold']), node)
            b_n = node
    d_out[b_n] = set()

    curr_path = [a_n]
    curr_v = a_n
    circuit = []

    while len(curr_path) > 0:
        if len(d_out[curr_v]) > 0:
            curr_path.append(curr_v)
            next_v = d_out[curr_v].pop()
            curr_v = next_v
        else:
            circuit.append(curr_v)
            curr_v = curr_path[-1]
            curr_path.pop()

    return circuit[::-1]


def createRandomSeq(L):
    seq = ""
    idx_to_nucleotide = {0:"A", 1:"T", 2:"C", 3:"G"}
    random.seed(datetime.now())
    for i in range(L):
        seq += idx_to_nucleotide[random.randint(0,3)]
    return seq

def seqToFasta(seq, seqName, filename):
    file = open(filename, "w")
    file.write(">"+seqName+"\n")
    file.write(seq)
    file.close()

def createRandomReads(seq, read_length, num_reads):
    reads = []
    idx_covered = set()
    max_idx = len(seq) - read_length - 1
    for i in range(num_reads):
        start = random.randint(0,max_idx)
        reads.append(seq[start:start+read_length+1])
        idx_covered.update(set(range(start,start+read_length+1)))

    # bin_length = int(read_length/2)
    # num_bins = int(max_idx/bin_length)
    # for i in range(num_bins):
    #     for _ in range(num_reads//num_bins):
    #         bin_start = bin_length*(i-1)
    #         bin_end = bin_length*i
    #         start = random.randint(bin_start,bin_end)
    #         reads.append(seq[start:start+read_length+1])
    #         idx_covered.update(set(range(start,start+read_length+1)))

    coverage = len(idx_covered)/len(seq)
    return reads, coverage

def createRandomSeq(L):
    seq = ""
    idx_to_nucleotide = {0:"A", 1:"T", 2:"C", 3:"G"}
    random.seed(datetime.now())
    for i in range(L):
        seq += idx_to_nucleotide[random.randint(0,3)]
    return seq

def merge(L, R):
    overlap = get_overlap(L,R)
    L += overlap
    return L

def get_overlap(s1, s2):
    m = min(len(s1), len(s2))
    for i in range(m,0,-1):
        if s1[-i:] == s2[:i]:
            return s2[i:]

def walk_to_sequence(walk):
    final = walk[0]
    for i in range(1,len(walk)):
        final = merge(final, walk[i])
    return final

def create_graph(reads):
    nodes = set()
    out_graph = dict()
    out_count = dict()
    in_graph = dict()
    in_count = dict()

    for i in range(len(reads)):
        L = reads[i][0:len(reads[i])-1]
        R = reads[i][1:]
        nodes.add(L)
        nodes.add(R)

        if L not in out_graph:
            out_graph[L] = {R}
        out_graph[L].add(R)

        if R not in in_graph:
            in_graph[R] = {L}
        in_graph[R].add(L)
        # print(L,"->",R)

    for n in nodes:
        if n not in out_graph:
            out_count[n] = 0
        else:
            out_count[n] = len(out_graph[n])

        if n not in in_graph:
            in_count[n] = 0
        else:
            in_count[n] = len(in_graph[n])

    return out_graph, out_count, in_graph, in_count


if __name__ == "__main__":

    print(colored("Generating Sequence", 'cyan', attrs=['bold']))
    seq = createRandomSeq(1000)

    print(colored("Generating Random Reads", 'cyan', attrs=['bold']))
    reads, coverage = createRandomReads(seq, 25, 100000)

    print(colored("Generating Graph", 'cyan', attrs=['bold']))
    out_graph, out_count, in_graph, in_count = create_graph(reads)

    print(colored("Performing Eulerian Walk", 'cyan', attrs=['bold']))
    start = timeit.default_timer()
    w  = perform_walk(in_graph, out_graph, in_count, out_count)
    stop = timeit.default_timer()
    print('Time: ', stop - start)

    print(colored("Generating Final Sequence", 'cyan', attrs=['bold']))
    final = walk_to_sequence(w)

    print(colored("seq   :", 'magenta', attrs=['bold']), seq)
    print(colored("result:", 'magenta', attrs=['bold']), final)

    if seq == final:
        print(colored("100% Match:", 'cyan', attrs=['bold']), colored('Yes', 'green', attrs=['bold']))
    else:
        print(colored("100% Match:", 'cyan', attrs=['bold']), colored('No', 'red', attrs=['bold']))

    print("Sequence Length and Runtime Experiment")
    runtimes = list()
    for i in tqdm(range(1000, 100000, 5000)):
        seq = createRandomSeq(i)
        reads, coverage = createRandomReads(seq, 25, 2000000)
        start = timeit.default_timer()
        og, oc, ig, ic = create_graph(reads)
        w  = perform_walk(ig, og, ic, oc)
        stop = timeit.default_timer()
        final = walk_to_sequence(w)
        runtimes.append(stop-start)
        if seq != final:
            print("fail")
    plt.figure(1)
    plt.plot(list(range(1000, 100000, 5000)), runtimes)
    plt.ylabel('Runtime (s)')
    plt.xlabel('Sequence Length (bp)')
    plt.savefig("runtime_seqlen.png", pad_inches=0.1)
    # plt.show()

    print("Runtime as Read Length is Varied")
    runtimes = list()
    for i in tqdm(range(20,1020,100)):
        seq = createRandomSeq(10000)
        reads, coverage = createRandomReads(seq, i, 200000)
        start = timeit.default_timer()
        og, oc, ig, ic = create_graph(reads)
        w  = perform_walk(ig, og, ic, oc)
        final = walk_to_sequence(w)
        stop = timeit.default_timer()
        runtimes.append(stop-start)
        if seq != final:
            print("fail")
    plt.figure(2)
    plt.plot(list(range(20,1020,100)), runtimes)
    plt.ylabel('Runtime (s)')
    plt.xlabel('Read Length (bp)')
    plt.savefig("runtime_readlength.png", pad_inches=0.1)
    # plt.show()

    print("Runtime as Number of Reads is Varied")
    runtimes = list()
    for i in range(10000, 400000, 5000):
        seq = createRandomSeq(1000)
        reads, coverage = createRandomReads(seq, 25, i)
        start = timeit.default_timer()
        og, oc, ig, ic = create_graph(reads)
        w  = perform_walk(ig, og, ic, oc)
        final = walk_to_sequence(w)
        stop = timeit.default_timer()
        runtimes.append(stop-start)
        if seq != final:
            print("fail")
    plt.figure(3)
    plt.plot(list(range(10000, 400000, 5000)), runtimes)
    plt.ylabel('Runtime (s)')
    plt.xlabel('Number of Reads')
    plt.savefig("runtime_numreads.png", pad_inches=0.1)
    # plt.show()

    print("Success Rate of Sequence Reconstruction based on Read Length")
    success_rates = []
    for i in tqdm(range(6,25,1)):
        successes = 0
        iter = 1000
        for j in range(iter):
            seq = createRandomSeq(1000)
            reads, coverage = createRandomReads(seq, i, 10000)
            og, oc, ig, ic = create_graph(reads)
            try:
                w  = perform_walk(ig, og, ic, oc)
                final = walk_to_sequence(w)
                if seq == final:
                    successes += 1
            except:
                continue
        success_rates.append(successes/iter)
        # print(successes/iter)
    plt.figure(4)
    plt.plot(list(range(6,25,1)), success_rates)
    plt.ylabel('Success Rate')
    plt.xlabel('Read Length (bp)')
    plt.savefig("successrate_readslen.png", pad_inches=0.1)
    plt.show()
