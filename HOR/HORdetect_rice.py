import networkx as nx
from networkx.drawing.nx_agraph import write_dot
import math
import os
from subprocess import check_call
import sys
import csv
import itertools
import argparse

def parse_args():
    parser = argparse.ArgumentParser(description='detect HOR on compressed dimer-seq')
    parser.add_argument('-fi', dest="fasta",help="dimer-compressed sequence [fasta-format]",required=False)
    parser.add_argument('-pos', dest="fpos",help="dimer genome coordinate [bed-format]",required=False)
    parser.add_argument('-p', dest="prefix",help="target sequence id for detection",required=False)
    parser.add_argument('-O', dest="outdir",help="directory for output file", default=".",required=False)
    parser.add_argument('-nodeThr', dest="nodeThr",help="remove monomers in graph with occurrences below the node threshold [default:2]", default=2,required=False)
    parser.add_argument('-edgeThr', dest="edgeThr",help="remove edges in graph with occurrence below the edge threshold [default : 2] ", default=2,required=False)
    parser.add_argument('-minTrav', dest="minTrav",help="minimum HOR occurance [default : 2] ", default=2,required=False)

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit()
    else:
        return parser.parse_args()


def stat_seq(fastafile):
    chrseq = {}
    with open(fastafile, 'r') as inf:
        name, seq = None, []
        for line in inf:
            line = line.rstrip()
            if line.startswith(">"):
                if name:
                    # yield (name[1:], ''.join(seq))
                    chrseq[name[1:]] = ''.join(seq)
                name, seq = line, []
            else:
                seq.append(line)
        if name:
            # yield (name[1:], ''.join(seq))
            chrseq[name[1:]] = ''.join(seq)

    return chrseq


def slidingWindow(sequence, winSize, step=1):

    out = []

    numOfChunks = ((len(sequence) - winSize) // step) + 1

    for i in range(0, numOfChunks * step, step):
        # yield i//step, sequence[i:i + winSize]
        out.append(sequence[i:i + winSize])

    return out


def stat_compression (ichr_mngroup, maxk, outprefix):
    for ichr, mngroup in ichr_mngroup.items():

        if ichr == outprefix:
        
            mngroups = slidingWindow(mngroup, 1, 1)
            ichr_k_cnt = [{} for k in range(maxk + 1)]
            
            for i, g in enumerate(mngroups):
                
                cur_mn = (g,)
                if cur_mn not in ichr_k_cnt[0]:
                    ichr_k_cnt[0][cur_mn] = 0
                ichr_k_cnt[0][cur_mn] += 1
                
                
                for k in range(1, maxk + 1):
                    
                    if i-k >= 0:
                        
                        cur_mn = (mngroups[i-k], *cur_mn)

                        if len(list(cur_mn)) == k+1:
                            if cur_mn not in ichr_k_cnt[k]:
                                ichr_k_cnt[k][cur_mn] = 0
                            ichr_k_cnt[k][cur_mn] +=1
                    else:
                        break
            return ichr_k_cnt


def BuildMonomerGraph(ichr_k_cnt, ichr, nodeThr=2, edgeThr=2 ):
    kcnt1 = ichr_k_cnt[0]
    kcnt2 = ichr_k_cnt[1]

    #print(kcnt2.keys())
    mn_set = {tuple(list(x)[:-1]) for x, y in kcnt2.items() if y > nodeThr} | \
        {tuple(list(x)[1:]) for x, y in kcnt2.items() if y > nodeThr}

    G = nx.DiGraph(label=ichr)

    for vert in mn_set:
        cnt_mn = 0

        if vert in kcnt1:
            cnt_mn = kcnt1[vert]
            lg = 0
            if cnt_mn > 0:
                lg = math.log(cnt_mn)

        clr = ["red", "#cd5700", "orange", "#ffb300", "yellow", "#ceff1d", "#a7fc00", "#00ff00", "#e0ffff", "#f5fffa",
           "#f5fffa", "#f5fffa", "#f5fffa", "#f5fffa"]
        vert = "-".join(list(vert))
        curc = clr[int(lg)]

        G.add_node(vert, style="filled", fillcolor=curc, label=r"%s\n[%s]" % (vert, str(int(cnt_mn))))

    ecnt = 0

    for vt1 in sorted(mn_set):
        for vt2 in sorted(mn_set):
            if list(vt1)[1:] != list(vt2)[:-1]:
                continue
            scr = 0
            if (*vt1, vt2[-1]) in kcnt2:
                scr = kcnt2[(*vt1, vt2[-1])]
            thr_wg = [100000000, 1000, 500, 100, 1]
            wgs = [7, 5, 3, 1, 0]
            wg = 3
            while (scr > thr_wg[wg] and wg > 0):
                 wg -= 1

            if scr > edgeThr:
                ecnt += 1
                vrt1 = "-".join(list(vt1))
                vrt2 = "-".join(list(vt2))
                if scr < 1:
                    G.add_edge(vrt1, vrt2, label=str(int(scr)), penwidth=str(wgs[wg]), constraint="false")
                else:
                    G.add_edge(vrt1, vrt2, label=str(int(scr)), penwidth=str(wgs[wg]))
    return G


def drawGraph(G, outprefix, outdir):
    write_dot(G, os.path.join(outdir, outprefix + ".dot"))
    check_call(['dot', '-Tpng', os.path.join(outdir, outprefix + ".dot"), '-o', os.path.join(outdir, outprefix + ".png")])


def addNode(G, mn, mncnt):
    lg = 0.01
    if mncnt > 0:
        lg = math.log(mncnt)
    clr = ["red", "#cd5700", "orange", "#ffb300", "yellow", "#ceff1d", "#a7fc00", "#00ff00", "#e0ffff", "#f5fffa"]
    G.add_node(mn, style=f'filled', fillcolor=f'{clr[int(lg)]}', label=f'{mn}[{str(int(mncnt))}]')


def addEdges(G, mn1, mn2, cnt, thr):
    scr = cnt
    thr_wg = [100000000, 1000, 500, 100, 1]
    wgs = [7, 5, 3, 1, 0]
    wg = 3
    while (scr > thr_wg[wg]):
        wg -= 1

    if scr > thr:
        G.add_edge(mn1, mn2, label=f'{scr}', penwidth=f'{wgs[wg]}')


def hor_non_continuity(cl, loci_index):
    horstr_len = len(cl) - 1
    if len(loci_index) == 0:
        return True
    else:
        dis = [loci_index[i+1] - loci_index[i] for i in range(len(loci_index)-1)]
        tflag = ""
        ### n(HORs) > 3 && distance(HOR1, HOR2) < len(HOR)*4
        for d in dis:
            if d < horstr_len * 4:
                tflag += "a"
            else:
                tflag += "c"
        # print(cl, ":", tflag)
        if "aa" in tflag:
            return False
        else:
            return True

def genCycleflag(cl, seq):
    clstr = "".join(cl)[:-1] # input cycle without the head to tail
    clstr_rotate = [clstr[i:] + clstr[:i] for i in range(len(clstr))]
    for iclstr in clstr_rotate:
        if seq.count(iclstr) >= 2:
            return  False
    else:
        return True

def genCycleInner(G, prefixCycle, usedEdges, cycleList, seq):
    def samecc(c1, c2):
        if len(c1) != len(c2):
            return False
        for j in range(len(c1) - 1):
            if c1[j:-1] + c1[:j] == c2[:-1]:
                return True
        return False

    if len(prefixCycle) > 1 and prefixCycle[0] == prefixCycle[-1]:
        # print(prefixCycle)
        for cc in cycleList:
            if samecc(cc, prefixCycle):
                break
            elif genCycleflag(prefixCycle, seq):
                break
        else:
            print("append cycles &: ",prefixCycle )
            cycleList.append(prefixCycle)
    v = prefixCycle[-1]
    for e in G.edges(v):
        if e not in usedEdges:
            flag = genCycleInner(G, prefixCycle + [e[1]], usedEdges | {e}, cycleList, seq)
            if flag:
                return flag
    return False


def genAllCycles(G, seq):
    cycleList = []
    for v in G.nodes():
        if G.in_degree(v) >= 1 and G.out_degree(v) >= 1:
            print("currect node: ", v)
            genCycleInner(G, [v], set(), cycleList, seq)
            # genCycleInner(G, [v], set(), cycleList)
    return cycleList


def extendcl(cl, maxclcnt ,maxindex, seq):
    out_clstr=""
    clstr = "".join(cl)[:-1]
    clstr_rotate, clstr_cnt, index = get_ClstrCnt(clstr, seq)
    extend_flag = False
    extend_rst = []
    for i in range(len(clstr_rotate)):
        for j in index[i]:
            while True:
                start = j+len(clstr)
                end = j+len(clstr)+1
                ss = seq[start: start + 2]
                if clstr[-2] == seq[start:end] and "".join(clstr_rotate[i])[:2] != ss:
                    update_clstr = seq[j:end]
                    update_cl = slidingWindow(update_clstr+update_clstr[0],1,1)
                    update_maxcl, update_maxclcnt, update_maxindex = getCycleCnt(update_cl, seq)
                    update_d = [update_maxindex[k + 1] - update_maxindex[k] for k in range(len(update_maxindex) - 1)]
                    cl_d = [maxindex[k + 1] - maxindex[k] for k in range(len(maxindex) - 1)]
                    if not hor_non_continuity(update_maxcl, update_maxindex):
                        if update_d.count(len(update_cl)-1) >= cl_d.count(len(cl)-1) or update_maxclcnt >= maxclcnt:
                            if update_maxindex != maxindex:
                                clstr = update_clstr
                                maxclcnt = update_maxclcnt
                                out_clstr = update_clstr
                                extend_flag = True
                                extend_rst.append(out_clstr)
                            else:
                                break
                        else:
                            break
                    else:
                        break
                else:
                    break

    if not extend_flag:
        return cl, maxclcnt, maxindex
    else:
        outcl = slidingWindow(out_clstr+out_clstr[0],1,1)
        update_maxcl, update_maxclcnt, update_maxindex = getCycleCnt(outcl, seq)
        update_d = [update_maxindex[k+1] - update_maxindex[k] for k in range(len(update_maxindex)-1)]
        cl_d = [maxindex[k+1] - maxindex[k] for k in range(len(maxindex)-1)]
        if cl_d.count(len(cl)-1) < update_d.count(len(update_maxcl)-1):
            return cl, maxclcnt, maxindex
        else:
            return update_maxcl, update_maxclcnt, update_maxindex

def extendCycle(cycles, cycles_cnt,cycles_index, seq):
    if len(cycles) != len(cycles_cnt) != len(cycles_index):
        print("something wrong for HOR index.")
        print(cycles,cycles_cnt,cycles_cnt )
        return
    outcycles = []
    outcycles_cnt = []
    outcycles_index = []
    for i in range(len(cycles)):
        cl = cycles[i]
        clcnt = cycles_cnt[i]
        clindex = cycles_index[i]
        maxcl, maxclcnt, maxindex = getCycleCnt(cl, seq)
        update_maxcl, update_maxclcnt, update_maxindex = extendcl(cl, maxclcnt, maxindex, seq)
        if update_maxcl == cl:
            outcycles.append(cl)
            outcycles_cnt.append(clcnt)
            outcycles_index.append(clindex)
        else:
            # outcycles.append([update_maxcl, cl])
            # outcycles_cnt.append([update_maxclcnt,clcnt])
            # outcycles_index.append([update_maxindex, clindex])
            outcycles.append(update_maxcl)
            outcycles_cnt.append(update_maxclcnt)
            outcycles_index.append(update_maxindex)
    return outcycles, outcycles_cnt, outcycles_index

def get_ClstrCnt(clstr, seq):
    # rotate the cycle string and get the rotate/cnt/index list #
    clstr_rotate = [clstr[i:] + clstr[:i] for i in range(len(clstr))]
    clstr_rotate = sorted(clstr_rotate)
    clstr_cnt = []
    index = []
    for iclstr in clstr_rotate:
        i = 0
        clcnt = 0
        clindex=[]
        while True:
            if i > len(seq):
                break
            if seq[i:i+len(iclstr)] == iclstr:
                clcnt += 1
                clindex.append(i)
                i = i+len(iclstr)
            else:
                i += 1
        # for i in range(len(seq) - len(iclstr)):
        #     if seq[i:i+len(iclstr)] == iclstr:
        #         clcnt += 1
        #         clindex.append(i)
        clstr_cnt.append(clcnt)
        index.append(clindex)
    return clstr_rotate, clstr_cnt, index

def getCycleCnt(cl, seq):
    # get the max Cnt after rotate the cycle with one base step #
    # input cl with head to tail #
    # returan cycle with head to tail #
    clstr = "".join(cl)[:-1]
    clstr_rotate, clstr_cnt, index = get_ClstrCnt(clstr, seq)
    maxclcnt = max(clstr_cnt)
    maxclstr = clstr_rotate[clstr_cnt.index(maxclcnt)]
    maxcl = slidingWindow(maxclstr+maxclstr[0], 1, 1)
    maxindex = index[clstr_cnt.index(maxclcnt)]
    return maxcl, maxclcnt, maxindex


def filterCycles(cycles, hybridSet, minTrav, seq):
    # return cl with head to tail #
    usedcycle = []
    cycles.sort(key=lambda x: -len(x))
    res_cyc = []
    res_cnt = []
    res_index = []

    def check_subHOR(cl, usecl):
        subhor_flag = False
        clstr_rotate, clstr_cnt, index = get_ClstrCnt("".join(cl)[:-1], seq)
        for rclstr in clstr_rotate:
            if rclstr in "".join(usecl)[:-1]*2:
                subhor_flag = True
        return subhor_flag

    for cl in cycles:
        if len([v for v in cl if v in hybridSet]) > 0:
            continue

        cl, clstr_cnt, clindex = getCycleCnt(cl,seq)
        # return cl with head to tail #
        if len(cl) <= 3:
            print("filter cycle:", cl, " for minNodes.")
            continue

        if clstr_cnt < minTrav:
            print("filter cycle: ", cl, " for minTrav, Trav is ", clstr_cnt)
            continue

        if hor_non_continuity(cl, clindex):
            print("filter cycle", cl, "num:", clstr_cnt, "index: " ,clindex )
            print("filter cycle: ", cl, " for non_continuity")
            continue

        filterflag = []
        removecycle = []
        if len(usedcycle) != 0:
            for usedc in usedcycle:
                if check_subHOR(cl, usedc):
                    usecl, use_cnt, use_index = getCycleCnt(usedc, seq)
                    use_d = [use_index[k+1] - use_index[k] for k in range(len(use_index)-1)]
                    cl_d = [clindex[k+1] - clindex[k] for k in range(len(clindex)-1)]
                    # iter_flag = any(sum(1 for _ in g) > 1 for _, g in itertools.groupby(cl_d) if _ == len(cl)-1)
                    if (cl_d.count(len(cl)-1) <= use_d.count(len(usecl)-1) + 1) or use_d.count(len(usecl)-1) >= 2 :
                        print("filter cycle", cl, "num:", clstr_cnt, "index: ", clindex)
                        print("filter cycle: ", cl, " for subHOR of ", usedc)
                        filterflag.append(True)
                    else:
                        print("cycle:", cl, "num", clstr_cnt, "index: ", clindex)
                        print("usecl:", usecl, "usenum", use_cnt, "index: ", use_index)
                        print(cl_d.count(len(cl)-1))
                        print(use_d.count(len(usecl)-1))
                        del res_cnt[res_cyc.index(usecl)]
                        del res_index[res_cyc.index(usecl)]
                        res_cyc.remove(usecl)
                        removecycle.append(usecl)
                        print("usecycles:", usedcycle)
            for recl in removecycle:
                usedcycle.remove(recl)
            if filterflag == []:
                usedcycle.append(cl)
                res_cyc.append(cl)
                res_cnt.append(clstr_cnt)
                res_index.append(clindex)
        else:
            usedcycle.append(cl)
            res_cyc.append(cl)
            res_cnt.append(clstr_cnt)
            res_index.append(clindex)
    return res_cyc, res_cnt, res_index


def merge_nestHOR(cycles, cycles_index, seq):
    out_cycles = []
    out_cycles_cnt = []
    out_cycles_index = []
    usedcycles = []
    tmp_cycles = []
    tmp_cycles_cnt = []
    tmp_cycles_index = []
    import numpy as np
    max_index = np.max([np.max(i) for i in cycles_index])
    max_len = np.max([len(c)-1 for c in cycles])

    max_count = np.max([len(i) for i in cycles_index])

    print(max_len, max_index)
    array = np.zeros(((max_index + max_len), ), dtype=int)
    num = 1
    for i in range(len(cycles)):
        delta = max_count
        for index in cycles_index[i]:
            for j in range(len(cycles[i])-1):
                if array[index + j] == 0:
                    array[index + j] = num
                else:
                    array[index + j] = -1
            num += 1
            delta -= 1
            # array[index: index + len(cycles[i])] = 1
        num += delta
    print(array)
    begin = 0
    end = 0
    start = False
    stop = False
    prev_num = 0
    rst = []
    check_is_same = np.zeros(len(cycles), dtype=int)
    for i in range(array.shape[0]):
        if array[i] == 0:
            if start:
                if stop:
                    start = False
                    stop = False
                else:
                    end = i
                    start = False
                    stop = False
                    rst.append([begin, end])
        else:
            if not start:
                check_is_same = np.zeros(len(cycles))
                start = True
                stop = False
                begin = i
                prev_num = array[i]
            if start and not stop:
                if array[i] == prev_num or prev_num == -1:
                    if prev_num == -1 and array[i] != prev_num:
                        if check_is_same[array[i] // max_count] == 1:
                            stop = True
                            end = i - 1
                            # start = False
                            rst.append([begin, end])
                    prev_num = array[i]
                    continue
                elif array[i] == -1:
                    if prev_num > 0:
                        print(prev_num, max_count)
                        check_is_same[prev_num // max_count] = 1
                    prev_num = -1
                    continue
                else:
                    start = False
                    end = i
                    rst.append([begin, end])

    if start:
        end = array.shape[0]
        rst.append([begin, end])

    # print(rst)
    for begin, end in rst:
        candidate_clstr = seq[begin: end]
        candidate_cl = slidingWindow(candidate_clstr+candidate_clstr[0],1,1)
        candidate_cl, candidate_clcnt, candidate_index = getCycleCnt(candidate_cl, seq)
        print("candidateclstr:", candidate_clstr)
        print("candidatecycle:", candidate_cl)
        if candidate_clcnt > 2 and not hor_non_continuity(candidate_cl, candidate_index) and candidate_cl not in tmp_cycles:
            tmp_cycles.append(candidate_cl)
            tmp_cycles_cnt.append(candidate_clcnt)
            tmp_cycles_index.append(candidate_index)
    combined = list(zip(tmp_cycles, tmp_cycles_cnt, tmp_cycles_index))
    sorted_combined = sorted(combined, key=lambda x: -len(x[0]))
    #sorted_cycles = sorted(tmp_cycles, key=lambda x: -len(x))
    sorted_cycles, sorted_cycles_cnt, sorted_cycles_index = zip(*sorted_combined)
    print("sorted_cycles:", sorted_cycles)
    print("sorted_cycles_cnt", sorted_cycles_cnt)
    print("sorted_cycles_index", sorted_cycles_index)
    for n, cl in enumerate(sorted_cycles):
        filterflag = False
        if usedcycles != []:
            for usedclstr in usedcycles:
                if "".join(cl)[:-1] in usedclstr:
                    cl, clcnt, index = getCycleCnt(cl, seq)
                    usedcl, usedcnt, usedindex = getCycleCnt(slidingWindow(usedclstr+usedclstr[0],1,1), seq)
                    use_d = [usedindex[k+1] - usedindex[k] for k in range(len(usedindex)-1)]
                    cl_d = [index[k+1] - index[k] for k in range(len(index)-1)]
                    # iter_flag = any(sum(1 for _ in g) > 1 for _, g in itertools.groupby(cl_d) if _ == len(cl)-1)
                    if (cl_d.count(len(cl)-1) <= use_d.count(len(usedcl)-1) + 1) or use_d.count(len(usedcl)-1) >= 2 :
                        filterflag = True
                        break
                    else:

                        _i = out_cycles.index(usedcl)
                        out_cycles_cnt.pop(_i)
                        out_cycles_index.pop(_i)
                        out_cycles.pop(_i)
            if not filterflag:
                out_cycles.append(cl)
                _i = sorted_cycles.index(cl)
                out_cycles_cnt.append(sorted_cycles_cnt[_i])
                out_cycles_index.append(sorted_cycles_index[_i])
                usedcycles.append("".join(cl)[:-1])
        else:
            out_cycles.append(cl)
            out_cycles_cnt.append(sorted_cycles_cnt[sorted_cycles.index(cl)])
            out_cycles_index.append(sorted_cycles_index[sorted_cycles.index(cl)])
            usedcycles.append("".join(cl)[:-1])

    return out_cycles, out_cycles_cnt, out_cycles_index

def detectHORs(G, minTrav, seq):
    hybridSet = {}
    # head to tail cycles
    cycles = genAllCycles(G, seq)
    print("AllCycles:", cycles)
    # print("AllCycles:", len(cycles))
    # head to tail cycles
    cycles, cycles_cnt, cycles_index = filterCycles(cycles, hybridSet, minTrav, seq)
    print("AfterFiterCycles:", cycles)
    print("AfterFiterCycles:", cycles_cnt)
    print("AfterFiterCycles:", cycles_index)
    cycles, cycles_cnt, cycles_index = extendCycle(cycles, cycles_cnt, cycles_index, seq)
    print("ExtendCycles:", cycles)
    print("ExtendCycles:", cycles_cnt)
    print("ExtendCycles:", cycles_index)
    if len(cycles) >1:
        cycles, cycles_cnt, cycles_index = merge_nestHOR(cycles, cycles_index, seq)
        print("MergeCycles:", cycles)
        print("MergeCycles index:", cycles_index)
    return cycles, cycles_cnt, cycles_index

def output_loci_index(horstr, loci_index):
    ### filter distance(HOR1, HOR2) > len(HOR)*4
    out = []
    tmp=[]
    horlen = len(horstr)
    # dis = [loci_index[i+1] - loci_index[i] for i in range(len(loci_index)-1)]
    for i in range(len(loci_index)-1):
        d = loci_index[i+1] - loci_index[i]
        if tmp == []:
            if  d < horlen * 4:
                tmp.append(loci_index[i])
                tmp.append(loci_index[i+1])
        else:
            if d < horlen * 4 and i != len(loci_index)-2:
                tmp.append(loci_index[i+1])
            elif d < horlen * 4 and i == len(loci_index)-2:
                tmp.append(loci_index[i + 1])
                out.append(tmp)
            else:
                out.append(tmp)
                tmp = []
    print(horstr, ":", out)
    inum = 0
    outindex = []
    flag = False
    for i, iindex in enumerate(out):
        if len(iindex) >= 3:
            flag = True
            inum += len(iindex)
            tmp = ",".join([str(j) for j in iindex])
            outindex.append(tmp)
    return flag, inum, outindex

def saveHOR(cycles, cycles_index, outprefix, outdir, ichr_dimer_pos):
    outfile = os.path.join(outdir,  outprefix+".HORs.tsv")
    with open(outfile, "w",  newline="") as fw:
        csv_writer = csv.writer(fw, delimiter='\t')
        if len(cycles) > 0 :
            for i in range(len(cycles)):
                cycle = cycles[i]
                horstr, loci_index = "".join(cycle)[:-1], cycles_index[i]
                horstr_rorate_sort = sorted([horstr[i:] + horstr[:i] for i in range(len(horstr))])
                ordered_horstr = horstr_rorate_sort[0]
                flag, inum, outindex = output_loci_index(horstr, loci_index)
                if flag:
                    for subindex in outindex:
                        pos_s, pos_e = stat_absolute_pos(horstr, subindex, ichr_dimer_pos, outprefix)
                        csv_writer.writerow(["H" + str(i+1), horstr, ordered_horstr, outprefix, len(subindex.split(',')), subindex, str(pos_s)+","+str(pos_e)])
        else:
            csv_writer.writerow(["-", "-", "-", outprefix, "-", "-"])
    return outfile

def decompose_sequence(seq, cycles, outprefix, outdir):
    outfile = os.path.join(outdir, outprefix + ".HORdecomposition.tsv")
    cycles = ["".join(cl)[:-1] for cl in cycles]
    with open(outfile, "w", newline="") as fw:
        csv_writer = csv.writer(fw, delimiter='\t')
        decomposed = []
        index = 0
        while index < len(seq):
            found_subtring = ""
            out_substring = ""
            terminate_flag = False
            found_substring_list = []
            cycle_substring_list = []
            single_base = []
            for cycle in cycles:
                cycle_rotate = [cycle[i:] + cycle[:i] for i in range(len(cycle))]
                if seq[index:index+len(cycle)] in cycle_rotate:
                    decomposed.append(seq[index:index+len(cycle)])
                    csv_writer.writerow([outprefix,  index, index+len(cycle), seq[index:index+len(cycle)]])
                    index += len(cycle)
                    break
                else:
                    i = 1
                    while True:
                        substring = seq[index:index + i]
                        if index + i == len(seq):
                            terminate_flag = True
                            break
                        if substring in cycle*2:
                            out_substring = substring
                            i += 1
                        elif len(out_substring) > 1:
                            found_substring_flag = True
                            found_substring_list.append(found_substring_flag)
                            cycle_substring_list.append(out_substring)
                            out_substring=""
                            break
                        else:
                            if len(substring) == 1:
                                single_base.append(substring)
                                break
                            else:
                                single_base.append(out_substring)
                                break
                    if found_substring_list != []:
                        for k in range(len(found_substring_list)):
                            if found_substring_list[k]:
                                if cycle_substring_list[k] > found_subtring:
                                    found_subtring = cycle_substring_list[k]
                        decomposed.append(found_subtring)
                        csv_writer.writerow([outprefix, index, index + len(found_subtring), found_subtring])
                        index += len(found_subtring)
                        break
                    elif terminate_flag:
                        decomposed.append(substring)
                        csv_writer.writerow([outprefix, index, index + len(substring), substring])
                        index += len(substring)
                        break
                    else:
                        if len(single_base) == len(cycles):
                            decomposed.append(single_base[0])
                            csv_writer.writerow([outprefix, index, index + len(single_base[0]), single_base[0]])
                            index += len(single_base[0])
                            break
    return outfile

def check_error(outprefix, outdir):
    flag = False
    infile = os.path.join(outdir, outprefix + ".HORs.tsv")
    with open(infile, 'r') as inf:
        for line in inf:
            tokens = line.strip().split("\t")
            if len(tokens) >1:
                horid, horstr, ordered_horstr, chrom, num, index, pos = tokens[:]
                if '-' not in index:
                    if ";" not in index:
                        if int(num) != len(index.split(",")):
                            flag = True
                            break
                    else:
                        tmp = [len(intervals.split(",")) for intervals in index.split(";")]
                        if int(num) != sum(tmp):
                            flag = True
                            break
        if flag:
            print("The index and numbers are inconsistent in ", outprefix)
        else:
            print("No errors in ", outprefix)

def read_absolute_pos(fpos):
    ichr_dimer_pos = {}
    with open(fpos, 'r') as inf:
        for line in inf:
            tokens = line.strip().split("\t")
            ichr = tokens[0]
            if ichr not in ichr_dimer_pos:
                ichr_dimer_pos[ichr] = [[int(tokens[1]), int(tokens[2])]]
            else:
                ichr_dimer_pos[ichr].append([int(tokens[1]), int(tokens[2])])
    return ichr_dimer_pos

def stat_absolute_pos(horstr, sub_index, ichr_dimer_pos, outprefix):
    dimer_pos = ichr_dimer_pos[outprefix]
    print(dimer_pos)
    cycle_sub_index = [int(i) for i in sub_index.split(",")]
    first_start = cycle_sub_index[0]
    last_end = cycle_sub_index[-1] + len(horstr)
    if first_start % 2 == 0:
        pos_s_index = first_start // 2
    else:
        pos_s_index = ( first_start + 1 ) // 2
    if last_end % 2 == 0:
        pos_e_index  = last_end // 2 - 1
    else:
        pos_e_index = last_end  // 2

    pos_s = dimer_pos[pos_s_index][0]
    pos_e = dimer_pos[pos_e_index][1]
    return pos_s, pos_e


def main():
    args = parse_args()
    fastafile = args.fasta
    posfile = args.fpos
    nodeThr = args.nodeThr
    edgeThr = args.edgeThr
    outprefix = args.prefix
    outdir = args.outdir
    minTrav = args.minTrav

    ichr_mngroup = stat_seq(fastafile)
    ichr_k_cnt = stat_compression (ichr_mngroup, 2, outprefix)
    G = BuildMonomerGraph(ichr_k_cnt, outprefix, nodeThr, edgeThr )
    drawGraph(G, outprefix, outdir)

    seq = ichr_mngroup[outprefix]
    print("Length of seq: ", len(seq))

    cycles, cycles_cnt, cycles_index = detectHORs(G, minTrav, ichr_mngroup[outprefix])
    ichr_dimer_pos = read_absolute_pos(posfile)
    outfile = saveHOR(cycles, cycles_index, outprefix, outdir, ichr_dimer_pos)
    print(outprefix, " detect HORs finished. saved in ", outfile)

    check_error(outprefix, outdir)
    
    if cycles !=[]:
        dec_outf = decompose_sequence(seq, cycles, outprefix, outdir)
        print(outprefix, " HOR decomposition finished. saved in ", dec_outf)



if __name__ == '__main__':
    main()
